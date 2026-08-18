# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This contains the logic for parsing MGI data files and producing Entry objects
for export to usable flat files.
"""

import csv
import logging
import typing as ty
from pathlib import Path

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.mgi.data import MgiEntry

LOGGER = logging.getLogger(__name__)

MOUSE_TAXID = 10090
REFSEQ_DBID = 9
ENSEMBL_DBID = 25

# MGI's 'Feature Type' is an SO term by another name, so this is a rename, not
# an inference. Anything mapping to None is not RNA and is dropped; anything
# absent is new since the last import and raises, because silently skipping a
# new ncRNA feature type is how a database quietly stops importing.
SO_TERMS: ty.Dict[str, ty.Optional[str]] = {
    "BAC/YAC end": None,
    "CpG island": None,
    "DNA segment": None,
    "QTL": None,
    "RNase MRP RNA gene": "SO:0000385",
    "RNase P RNA gene": "SO:0000386",
    "SRP RNA gene": "SO:0000590",
    "antisense lncRNA gene": "SO:0001904",
    "bidirectional promoter lncRNA gene": "SO:0001877",
    "chromosomal deletion": None,
    "complex/cluster/region": None,
    "endogenous retroviral region": None,
    "enhancer": None,
    "gene": None,
    "gene segment": None,
    "heritable phenotypic marker": None,
    "imprinting control region": None,
    "insertion": None,
    "intronic lncRNA gene": "SO:0001877",
    "lincRNA gene": "SO:0001877",
    "lncRNA gene": "SO:0001877",
    "locus control region": None,
    "miRNA gene": "SO:0000276",
    "minisatellite": None,
    "non-coding RNA gene": "SO:0000655",
    "open chromatin region": None,
    "origin of replication": None,
    "other genome feature": None,
    "polymorphic pseudogene": None,
    "promoter": None,
    "protein coding gene": None,
    "pseudogene": None,
    "pseudogenic gene segment": None,
    "pseudogenic region": None,
    "rRNA gene": "SO:0000252",
    "retrotransposon": None,
    "ribozyme gene": "SO:0000374",
    "scRNA gene": "SO:0000013",
    "sense intronic lncRNA gene": "SO:0001877",
    "sense overlapping lncRNA gene": "SO:0001877",
    "snRNA gene": "SO:0000274",
    "snoRNA gene": "SO:0000275",
    "tRNA gene": "SO:0000253",
    "telomerase RNA gene": "SO:0000390",
    "transcriptional cis regulatory region": None,
    "transgene": None,
    "unclassified cytogenetic marker": None,
    "unclassified gene": None,
    "unclassified non-coding RNA gene": "SO:0000655",
    "unclassified other genome feature": None,
}


def column_name(header: str) -> str:
    return header.strip().lower().replace(" ", "_")


def load(path: Path) -> ty.List[MgiEntry]:
    """
    Load every marker in an MRK_Sequence.rpt file.
    """

    with path.open("r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header in {path}")
        reader.fieldnames = [column_name(f) for f in reader.fieldnames]
        return [MgiEntry.from_row(row) for row in reader]


def so_term(entry: MgiEntry) -> ty.Optional[str]:
    """
    Get the SO term for a marker, or None if the marker is not RNA.
    """

    if entry.feature_type not in SO_TERMS:
        raise ValueError(f"Unknown MGI feature type: {entry.feature_type}")
    return SO_TERMS[entry.feature_type]


def url(entry: MgiEntry) -> str:
    return f"https://www.informatics.jax.org/marker/{entry.mgi_id}"


def description(entry: MgiEntry) -> str:
    return f"Mus musculus (house mouse) {entry.name}"


def references() -> ty.List[Reference]:
    return [
        Reference(
            authors=(
                "Blake JA, Eppig JT, Kadin JA, Richardson JE, Smith CL, Bult CJ; "
                "the Mouse Genome Database Group."
            ),
            location="Nucleic Acids Res. 2017 Jan 4;",
            title=(
                "Mouse Genome Database (MGD)-2017: community knowledge resource "
                "for the laboratory mouse"
            ),
            pmid=27899570,
            doi="10.1093/nar/gkw1040",
        )
    ]


def longest(rows) -> ty.Dict[str, str]:
    """
    Collapse (key, urs, length) rows to {key: urs}, keeping the longest
    sequence for each key.
    """

    best: ty.Dict[str, ty.Tuple[str, int]] = {}
    for key, urs, length in rows:
        if key not in best or length > best[key][1]:
            best[key] = (urs, length)
    return {key: urs for key, (urs, _) in best.items()}


def transcript_mapping(conn, transcript_ids: ty.List[str], dbid: int):
    """
    Map every transcript id to the longest mouse URS it corresponds to.

    Both RefSeq and Ensembl import the transcript id as external_id, but
    optional_id is checked too, as which of the two a database fills in has
    moved around. Scoping to one dbid matters: without it the planner probes
    all ~180 xref partitions, which turns seconds into minutes.
    """

    if not transcript_ids:
        return {}

    ids = sorted(set(transcript_ids))
    query = """
    select k, rna.urs, rna.len
    from (
        select acc.external_id as k, acc.accession
        from rnc_accessions acc where acc.external_id = ANY(%s)
      union all
        select acc.optional_id, acc.accession
        from rnc_accessions acc where acc.optional_id = ANY(%s)
    ) m
    join xref on xref.ac = m.accession
    join rna on rna.urs = xref.urs
    where xref.taxid = %s and xref.deleted = 'N' and xref.dbid = %s
    """
    with conn.cursor() as cur:
        cur.execute(query, (ids, ids, MOUSE_TAXID, dbid))
        return longest(cur)


def refseq_mapping(conn, refseq_ids: ty.List[str]) -> ty.Dict[str, str]:
    return transcript_mapping(conn, refseq_ids, REFSEQ_DBID)


def ensembl_mapping(conn, ensembl_ids: ty.List[str]) -> ty.Dict[str, str]:
    return transcript_mapping(conn, ensembl_ids, ENSEMBL_DBID)


def sequence_mapping(conn, urs_ids: ty.List[str]) -> ty.Dict[str, str]:
    """
    Fetch the sequence for every URS we resolved.
    """

    if not urs_ids:
        return {}

    query = """
    select urs, coalesce(seq_short, seq_long)
    from rna where urs = ANY(%s)
    """
    with conn.cursor() as cur:
        cur.execute(query, (sorted(set(urs_ids)),))
        return dict(cur)
