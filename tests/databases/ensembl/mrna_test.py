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

from pathlib import Path

import pytest

from rnacentral_pipeline.databases.data.regions import Strand
from rnacentral_pipeline.databases.ensembl import mrna

GFF = Path("data/ensembl/Homo_sapiens.GRCh38.mrna.gff3")


@pytest.fixture(scope="module")
def entries():
    with open("data/ensembl/Homo_sapiens.GRCh38.mrna.embl", "r") as raw:
        return {e.accession: e for e in mrna.parse(raw, GFF)}


def exons(entry):
    return [(e.start, e.stop) for e in entry.regions[0].exons]


def related(entry):
    return [(r.sequence_id, r.relationship) for r in entry.related_sequences]


def test_builds_one_entry_per_mrna_and_utr(entries):
    types = [e.rna_type for e in entries.values()]
    assert types.count("SO:0000234") == 17
    assert types.count("SO:0000204") == 4
    assert types.count("SO:0000205") == 4
    assert all(e.database == "ENSEMBL_MRNA" for e in entries.values())
    assert all(e.is_valid() for e in entries.values())


def test_builds_the_mrna(entries):
    entry = entries["ENST00000615165.1"]
    assert entry.primary_id == "ENST00000615165"
    assert entry.seq_version == "1"
    assert entry.ncbi_tax_id == 9606
    assert entry.feature_name == "mRNA"
    assert entry.description == "Homo sapiens (human) ENSG00000277196 mRNA"
    assert len(entry.sequence) == 1990
    assert entry.url == (
        "https://www.ensembl.org/feature-explorer/GCA_000001405.29/"
        "transcript:ENST00000615165"
    )


def test_spliced_utrs_on_the_minus_strand_match_ensembl(entries):
    five = entries["ENST00000615165.1:five_prime_UTR"]
    three = entries["ENST00000615165.1:three_prime_UTR"]
    assert five.regions[0].strand == Strand.reverse
    assert exons(five) == [(156447, 156497), (161689, 161750)]
    assert exons(three) == [(138082, 138479)]
    assert five.feature_name == "5'UTR"
    assert three.feature_name == "3'UTR"
    assert five.description == "Homo sapiens (human) ENSG00000277196 5' UTR"


def test_utrs_are_the_ends_of_the_mrna(entries):
    for accession in ["ENST00000615165.1", "ENST00000617983.1"]:
        mrna_seq = entries[accession].sequence
        assert mrna_seq.startswith(entries[f"{accession}:five_prime_UTR"].sequence)
        assert mrna_seq.endswith(entries[f"{accession}:three_prime_UTR"].sequence)


def test_links_the_mrna_and_utrs_both_ways(entries):
    assert related(entries["ENST00000617983.1"]) == [
        ("ENST00000617983.1:five_prime_UTR", "five_prime_utr"),
        ("ENST00000617983.1:three_prime_UTR", "three_prime_utr"),
    ]
    assert related(entries["ENST00000617983.1:five_prime_UTR"]) == [
        ("ENST00000617983.1", "mrna")
    ]
    assert related(entries["ENST00000617983.1:three_prime_UTR"]) == [
        ("ENST00000617983.1", "mrna")
    ]


def test_mitochondrial_mrnas_have_no_utrs(entries):
    entry = entries["ENST00000361390.2"]
    assert entry.description == "Homo sapiens (human) MT-ND1 mRNA"
    assert related(entry) == []
    assert "ENST00000361390.2:five_prime_UTR" not in entries


def test_only_protein_coding_transcripts_are_kept_with_their_tags(tmp_path):
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "1\thavana\tmRNA\t1\t9\t.\t+\t.\tID=transcript:ENST1;Parent=gene:ENSG1;"
        "biotype=protein_coding;tag=gencode_basic,MANE_Select,Ensembl_canonical\n"
        "1\thavana\tmRNA\t1\t9\t.\t+\t.\tID=transcript:ENST2;Parent=gene:ENSG1;"
        "biotype=nonsense_mediated_decay\n"
        "1\thavana\tlnc_RNA\t1\t9\t.\t+\t.\tID=transcript:ENST3;Parent=gene:ENSG2;"
        "biotype=protein_coding\n"
    )
    assert mrna.protein_coding(gff) == {
        "ENST1": ("ENSG1", ["Ensembl_canonical", "MANE_Select"])
    }


def test_the_mrna_keeps_its_protein_xrefs_and_tags(entries):
    entry = entries["ENST00000361390.2"]
    assert entry.xref_data == {
        "Uniprot/SWISSPROT": ["P03886"],
        "RefSeq_peptide": ["YP_003024026"],
        "Uniprot/SPTREMBL": ["U5Z754"],
        "HGNC": ["HGNC:7455"],
    }
    assert entry.note_data["tags"] == ["Ensembl_canonical"]
    assert entries["ENST00000615165.1"].note_data == {}


def test_utrs_carry_no_protein_xrefs(entries):
    utr = entries["ENST00000617983.1:five_prime_UTR"]
    assert utr.xref_data == {}
    assert utr.note_data == {}
