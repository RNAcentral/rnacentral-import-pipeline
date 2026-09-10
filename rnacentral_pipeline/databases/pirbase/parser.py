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

piRBase v3.0 ships one FASTA per species with nothing in the header but the
piRNA id, so the taxid, the rna_type and the description are all inferred here.

The full sets are far too large to import whole -- tens of millions of
sequences across the forty three species -- so they are filtered to piRNAs
whose sequence RNAcentral already holds. That keeps everything imported
previously, since those sequences are in `rna`, and admits new piRNAs only
where another database has deposited the same sequence independently.

The six gold standard sets are curated and small, so they bypass the filter
entirely. They have to: their ids barely overlap the full sets, and a new gold
sequence has nothing to match against.

No coordinates are imported. piRBase publishes BED against GRCm38, Rnor6 and
macFas5, none of which RNAcentral carries any more, and RNAcentral's own genome
mapping already places these sequences on the current assemblies.
"""

import logging
import typing as ty
from pathlib import Path

from Bio import SeqIO
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.pirbase.data import Species

LOGGER = logging.getLogger(__name__)

# Everything piRBase publishes is a piRNA; it has no other rna_type.
RNA_TYPE = "SO:0001035"


def build_known(md5s: Path, output: Path) -> None:
    """
    Index the md5s of every sequence already in RNAcentral. There are tens of
    millions of them, so this goes to sqlite rather than a set, and it is built
    once for the whole run rather than per species.

    The path is explicit because the parse tasks run under `--contain`, where
    the default tempfile location is a tmpfs of a few tens of MB.
    """

    known = SqliteDict(str(output))
    with md5s.open("r") as raw:
        for line in raw:
            known[line.strip()] = True
    known.commit()
    known.close()


def description(species: Species, accession: str) -> str:
    return f"{species.name} piRNA {accession}"


def as_entry(species: Species, record) -> data.Entry:
    return data.Entry(
        primary_id=record.id,
        accession=record.id,
        ncbi_tax_id=species.taxid,
        database="PIRBASE",
        sequence=str(record.seq).upper().replace("U", "T"),
        regions=[],
        rna_type=RNA_TYPE,
        # piRBase has no per-piRNA page; browse.php ignores a name parameter.
        url="",
        seq_version="1",
        description=description(species, record.id),
        species=species.name,
    )


def parse(
    handle: ty.IO, species: Species, known_path: ty.Optional[Path] = None
) -> ty.Iterable[data.Entry]:
    """
    Parse one piRBase FASTA. Without known_path every sequence is yielded, which
    is what the gold standard sets want; with it only sequences RNAcentral
    already holds survive.
    """

    known = SqliteDict(str(known_path), flag="r") if known_path else None

    total = 0
    kept = 0
    for record in SeqIO.parse(handle, "fasta"):
        total += 1
        entry = as_entry(species, record)
        if known is not None and entry.md5() not in known:
            continue
        kept += 1
        yield entry

    if known is not None:
        known.close()
    LOGGER.info("Parsed %i %s piRNAs, kept %i", total, species.code, kept)
