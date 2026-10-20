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

NONCODE ships one FASTA of lncRNA transcripts per species, with nothing in the
header but the versioned transcript id, and a BED12 of exon structure for the
two species it has genome coordinates for.
"""

import logging
import typing as ty
from pathlib import Path

from Bio import SeqIO

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.noncode.data import Species

LOGGER = logging.getLogger(__name__)

# Every sequence in these files is a lncRNA; NONCODE has no other rna_type.
RNA_TYPE = "SO:0001877"


def transcript_id(accession: str) -> str:
    """
    Strip the version, which is the id NONCODE's own pages are keyed on.
    """
    return accession.split(".")[0]


def seq_version(accession: str) -> str:
    parts = accession.split(".")
    return parts[1] if len(parts) > 1 else "1"


def url(accession: str) -> str:
    return f"http://www.noncode.org/show_rna.php?id={transcript_id(accession)}"


def description(species: Species, accession: str) -> str:
    return f"{species.name} long non-coding RNA {transcript_id(accession)}"


def regions(species: Species, bed_line: str) -> ty.List[data.SequenceRegion]:
    """
    Turn one BED12 line into a region. BED is half open and zero based, and the
    block starts are relative to the feature, not the chromosome.
    """

    fields = bed_line.split("\t")
    chromosome, start, _, _, _, strand = fields[:6]
    if chromosome.startswith("chr"):
        chromosome = chromosome[3:]
    sizes = [int(s) for s in fields[10].rstrip(",").split(",")]
    starts = [int(s) for s in fields[11].rstrip(",").split(",")]

    exons = [
        data.Exon(start=int(start) + offset, stop=int(start) + offset + size)
        for offset, size in zip(starts, sizes)
    ]

    return [
        data.SequenceRegion(
            assembly_id=species.assembly,
            chromosome=chromosome,
            strand=data.Strand.build(strand),
            exons=exons,
            coordinate_system=data.CoordinateSystem.zero_based(),
        )
    ]


def load_bed(handle: ty.Optional[ty.IO]) -> ty.Dict[str, str]:
    """
    Index a BED12 by the transcript id in its name column.
    """

    if handle is None:
        return {}
    return {line.split("\t")[3]: line.rstrip("\n") for line in handle if line.strip()}


def as_entry(
    species: Species, record, bed: ty.Dict[str, str]
) -> ty.Optional[data.Entry]:
    accession = record.id
    sequence = str(record.seq).upper().replace("U", "T")

    bed_line = bed.get(accession)
    return data.Entry(
        primary_id=transcript_id(accession),
        accession=accession,
        ncbi_tax_id=species.taxid,
        database="NONCODE",
        sequence=sequence,
        regions=regions(species, bed_line) if bed_line else [],
        rna_type=RNA_TYPE,
        url=url(accession),
        seq_version=seq_version(accession),
        description=description(species, accession),
        species=species.name,
    )


def parse(
    handle: ty.IO, species: Species, bed_handle: ty.Optional[ty.IO] = None
) -> ty.Iterable[data.Entry]:
    bed = load_bed(bed_handle)
    if bed and not species.assembly:
        raise ValueError(f"{species.key} has no assembly to place a BED against")

    placed = 0
    total = 0
    for record in SeqIO.parse(handle, "fasta"):
        entry = as_entry(species, record, bed)
        total += 1
        if entry.regions:
            placed += 1
        yield entry
    LOGGER.info("Parsed %i %s transcripts, %i placed", total, species.key, placed)
