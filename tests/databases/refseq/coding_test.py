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
"""

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from rnacentral_pipeline.databases.refseq import coding

DATA = "data/refseq/coding/{}.gbff"


def entries(name, taxid, assembly=None):
    with open(DATA.format(name)) as raw:
        return {e.accession: e for e in coding.parse(raw, taxid, assembly)}


def record(name):
    with open(DATA.format(name)) as raw:
        return next(SeqIO.parse(raw, "genbank"))


def test_a_virus_segment_is_one_mrna_with_utrs_around_all_its_orfs():
    found = entries("viruses", 11320)
    mrna = found["REFSEQ_CODING:NC_002020.1"]
    five = found["REFSEQ_CODING:NC_002020.1:five_prime_UTR"]
    three = found["REFSEQ_CODING:NC_002020.1:three_prime_UTR"]
    # NS1 is 27..719 and NEP join(27..56,529..864)
    assert mrna.sequence == str(record("viruses").seq)
    assert five.sequence == mrna.sequence[:26]
    assert three.sequence == mrna.sequence[864:]
    assert mrna.rna_type == "SO:0000234"
    assert five.rna_type == "SO:0000204"
    assert three.rna_type == "SO:0000205"
    assert mrna.description == (
        "Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 8 mRNA"
    )


def test_virus_entries_take_the_given_taxid_and_have_no_coordinates():
    found = entries("viruses", 11320)
    assert {e.ncbi_tax_id for e in found.values()} == {11320}
    assert all(not e.regions for e in found.values())


def test_virus_utrs_and_mrna_are_linked_both_ways():
    found = entries("viruses", 40051)
    mrna = found["REFSEQ_CODING:NC_006007.1"]
    assert {r.sequence_id for r in mrna.related_sequences} == {
        "REFSEQ_CODING:NC_006007.1:five_prime_UTR",
        "REFSEQ_CODING:NC_006007.1:three_prime_UTR",
    }
    utr = found["REFSEQ_CODING:NC_006007.1:three_prime_UTR"]
    assert [r.sequence_id for r in utr.related_sequences] == [mrna.accession]
    assert len(found["REFSEQ_CODING:NC_006007.1:five_prime_UTR"].sequence) == 19
    assert len(utr.sequence) == 41


def test_a_spliced_minus_strand_yeast_mrna_is_placed_on_the_assembly():
    found = entries("yeast", 28985, "ASM251v1")
    assert list(found) == ["REFSEQ_CODING:XM_002999299.1"]
    mrna = found["REFSEQ_CODING:XM_002999299.1"]
    genome = str(record("yeast").seq)
    assert mrna.sequence == str(
        Seq(genome[132:361] + genome[439:474]).reverse_complement()
    )
    assert mrna.sequence.startswith("ATG")
    [region] = mrna.regions
    assert region.assembly_id == "ASM251v1"
    assert region.chromosome == "A"
    assert region.strand.display_int() == -1
    assert [(e.start, e.stop) for e in region.exons] == [(133, 361), (440, 474)]
    assert mrna.gene == "PCC1"
    assert mrna.url == "https://www.ncbi.nlm.nih.gov/nuccore/XM_002999299.1"


def test_bacterial_cds_stand_in_for_mrnas_and_pseudogenes_are_skipped():
    found = entries("bacteria", 264730, "ASM1220v1")
    assert list(found) == ["REFSEQ_CODING:PSPPH_RS00030"]
    mrna = found["REFSEQ_CODING:PSPPH_RS00030"]
    assert mrna.regions[0].chromosome == "Chromosome"
    assert [(e.start, e.stop) for e in mrna.regions[0].exons] == [(73, 1398)]
    assert mrna.url == "https://www.ncbi.nlm.nih.gov/gene/?term=PSPPH_RS00030"


def test_organelle_records_are_skipped():
    # Ensembl has no mitochondrion for these yeasts, so there is nowhere to place them
    assert entries("mitochondrion", 28985, "ASM251v1") == {}


@pytest.mark.parametrize(
    "qualifiers,expected",
    [
        ({"chromosome": ["B"]}, "B"),
        ({"plasmid": ["large"]}, "large_plasmid"),
        ({}, "Chromosome"),
    ],
)
def test_chromosome_names_follow_ensembl(qualifiers, expected):
    rec = SeqRecord(Seq("ACGT"), id="NC_1.1")
    rec.features = [
        SeqFeature(FeatureLocation(0, 4), type="source", qualifiers=qualifiers)
    ]
    assert coding.chromosome(rec) == expected
