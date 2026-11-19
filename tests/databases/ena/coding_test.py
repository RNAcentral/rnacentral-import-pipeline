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

from rnacentral_pipeline.databases.ena import coding

DATA = "data/ena/coding/{}.embl"


def entries(name, taxid, assembly=None):
    with open(DATA.format(name)) as raw:
        return {e.accession: e for e in coding.parse(raw, taxid, assembly)}


def record(name):
    with open(DATA.format(name)) as raw:
        return next(SeqIO.parse(raw, "embl"))


def test_a_virus_segment_is_one_mrna_with_utrs_around_all_its_orfs():
    found = entries("viruses", 11320)
    mrna = found["ENA_CODING:J02150.1"]
    five = found["ENA_CODING:J02150.1:five_prime_UTR"]
    three = found["ENA_CODING:J02150.1:three_prime_UTR"]
    # NS1 is 27..719 and NEP join(27..56,529..864)
    assert mrna.sequence == str(record("viruses").seq)
    assert five.sequence == mrna.sequence[:26]
    assert three.sequence == mrna.sequence[864:]
    assert mrna.rna_type == "SO:0000234"
    assert five.rna_type == "SO:0000204"
    assert three.rna_type == "SO:0000205"
    assert mrna.description == (
        "Influenza A virus (A/Puerto Rico/8/34(H1N1)) segment 8 mRNA"
    )
    assert mrna.url == "https://www.ebi.ac.uk/ena/browser/view/J02150.1"


def test_virus_entries_take_the_given_taxid_and_have_no_coordinates():
    found = entries("viruses", 11320)
    assert {e.ncbi_tax_id for e in found.values()} == {11320}
    assert all(not e.regions for e in found.values())


def test_entries_cite_their_record_references():
    mrna = entries("viruses", 11320)["ENA_CODING:J02150.1"]
    assert "7465426" in {r.external_id for r in mrna.references}


# HIV's record is the whole pNL4-3 plasmid; only the provirus before the vector is viral
def test_a_record_can_be_cut_before_its_non_viral_end():
    with open(DATA.format("viruses")) as raw:
        found = {e.accession: e for e in coding.parse(raw, 11320, end=880)}
    mrna = found["ENA_CODING:J02150.1"]
    three = found["ENA_CODING:J02150.1:three_prime_UTR"]
    assert mrna.sequence == str(record("viruses").seq)[:880]
    assert three.sequence == mrna.sequence[864:]
    assert mrna.mol_type == "genomic RNA"
    assert "7465426" in {r.external_id for r in mrna.references}


def test_virus_utrs_and_mrna_are_linked_both_ways():
    found = entries("viruses", 40051)
    mrna = found["ENA_CODING:FJ969726.1"]
    assert {r.sequence_id for r in mrna.related_sequences} == {
        "ENA_CODING:FJ969726.1:five_prime_UTR",
        "ENA_CODING:FJ969726.1:three_prime_UTR",
    }
    utr = found["ENA_CODING:FJ969726.1:three_prime_UTR"]
    assert [r.sequence_id for r in utr.related_sequences] == [mrna.accession]
    assert len(found["ENA_CODING:FJ969726.1:five_prime_UTR"].sequence) == 19
    assert len(utr.sequence) == 41


def test_a_spliced_minus_strand_yeast_cds_is_placed_on_the_assembly():
    found = entries("yeast", 28985, "ASM251v1")
    assert list(found) == ["ENA_CODING:KLLA0_A00539g"]
    mrna = found["ENA_CODING:KLLA0_A00539g"]
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
    assert mrna.locus_tag == "KLLA0_A00539g"
    assert mrna.url == "https://www.ebi.ac.uk/ena/browser/view/CAR65216.1"


def test_bacterial_cds_stand_in_for_mrnas_and_pseudogenes_are_skipped():
    found = entries("bacteria", 264730, "ASM1220v1")
    assert list(found) == ["ENA_CODING:PSPPH_0005"]
    mrna = found["ENA_CODING:PSPPH_0005"]
    assert mrna.regions[0].chromosome == "Chromosome"
    assert [(e.start, e.stop) for e in mrna.regions[0].exons] == [(73, 1398)]
    assert mrna.url == "https://www.ebi.ac.uk/ena/browser/view/AAZ34480.1"


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
    rec = SeqRecord(Seq("ACGT"), id="CP000059.1")
    rec.features = [
        SeqFeature(FeatureLocation(0, 4), type="source", qualifiers=qualifiers)
    ]
    assert coding.chromosome(rec) == expected
