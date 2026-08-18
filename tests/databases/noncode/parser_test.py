# -*- coding: utf-8 -*-

"""
NONCODE gives us a bare versioned transcript id and a sequence, and for human
and mouse a BED12 of exon structure. Everything else -- taxid, rna_type, url,
description -- is inferred here, so this pins the inference.

The BED arithmetic is the part worth guarding: block starts are relative to the
feature start, not the chromosome, and BED is half open and zero based while
RNAcentral stores one based coordinates.
"""

from pathlib import Path

import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.noncode import parser
from rnacentral_pipeline.databases.noncode.data import SPECIES, species

YEAST = Path("data/noncode/yeast.fa")
HUMAN = Path("data/noncode/human.fa")
HUMAN_BED = Path("data/noncode/human.bed")


def parse(path, key, bed=None):
    with path.open("r") as handle:
        if bed is None:
            return list(parser.parse(handle, species(key)))
        with bed.open("r") as bed_handle:
            return list(parser.parse(handle, species(key), bed_handle))


def test_builds_entries_without_coordinates():
    entries = parse(YEAST, "yeast")
    assert len(entries) == 10

    entry = entries[0]
    assert entry.accession == "NONSCET000001.2"
    assert entry.primary_id == "NONSCET000001"
    assert entry.seq_version == "2"
    assert entry.ncbi_tax_id == 4932
    assert entry.database == "NONCODE"
    assert entry.rna_type == "SO:0001877"
    assert entry.url == "http://www.noncode.org/show_rna.php?id=NONSCET000001"
    assert entry.description == (
        "Saccharomyces cerevisiae long non-coding RNA NONSCET000001"
    )
    assert entry.regions == []
    assert entry.is_valid()


def test_places_transcripts_from_the_bed():
    entries = {e.accession: e for e in parse(HUMAN, "human", HUMAN_BED)}
    region = entries["NONHSAT000002.2"].regions[0]

    assert region.assembly_id == "GRCh38"
    # 'chr' is stripped, as everywhere else in the pipeline.
    assert region.chromosome == "1"
    assert region.strand == data.Strand.forward
    # chromStart 11871 with blockSizes 356,109,1188 at offsets 0,741,1353.
    assert [(e.start, e.stop) for e in region.exons] == [
        (11871, 12227),
        (12612, 12721),
        (13224, 14412),
    ]
    assert [(e.start, e.stop) for e in region.as_one_based().exons] == [
        (11872, 12227),
        (12613, 12721),
        (13225, 14412),
    ]


def test_every_transcript_in_the_bed_is_placed():
    entries = parse(HUMAN, "human", HUMAN_BED)
    assert all(e.regions for e in entries)


def test_transcripts_missing_from_the_bed_are_still_imported():
    entries = {e.accession: e for e in parse(HUMAN, "human", HUMAN_BED)}
    with HUMAN_BED.open("r") as handle:
        placed = set(parser.load_bed(handle))
    assert set(entries) >= placed


def test_refuses_a_bed_for_a_species_with_no_assembly():
    with HUMAN.open("r") as handle, HUMAN_BED.open("r") as bed:
        with pytest.raises(ValueError):
            list(parser.parse(handle, species("yeast"), bed))


def test_rejects_unknown_species():
    with pytest.raises(ValueError):
        species("narwhal")


def test_sequences_are_dna():
    entries = parse(YEAST, "yeast")
    assert all("U" not in e.sequence for e in entries)


@pytest.mark.parametrize("key", sorted(SPECIES))
def test_every_species_has_a_plausible_taxid(key):
    entry = SPECIES[key]
    assert entry.taxid > 0
    assert entry.name[0].isupper()
    # Only the two species with coordinates may carry an assembly.
    assert (entry.assembly is not None) == (key in {"human", "mouse"})
