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

from Bio import SeqIO
from Bio.Seq import Seq

from rnacentral_pipeline.databases.sgd import coding

GFF = "data/sgd/coding.gff"


def genome():
    with open(GFF) as raw:
        for line in raw:
            if line.startswith("##FASTA"):
                break
        return {r.id: str(r.seq) for r in SeqIO.parse(raw, "fasta")}


def piece(chrom, *spans, minus=False):
    seq = "".join(genome()[chrom][start - 1 : stop] for start, stop in spans)
    return str(Seq(seq).reverse_complement()) if minus else seq


def entries():
    with open(GFF) as raw:
        return {e.primary_id: e for e in coding.parse(raw)}


def test_builds_mrnas_and_utrs_from_the_mrna_span():
    found = entries()
    assert found["YAA001W_id001"].sequence == piece("chrI", (11, 80))
    assert found["YAA001W_id001:five_prime_UTR"].sequence == piece("chrI", (11, 20))
    assert found["YAA001W_id001:three_prime_UTR"].sequence == piece("chrI", (61, 80))
    assert found["YAA001W_id001"].rna_type == "SO:0000234"
    assert found["YAA001W_id001:five_prime_UTR"].rna_type == "SO:0000204"
    assert found["YAA001W_id001:three_prime_UTR"].rna_type == "SO:0000205"


def test_removes_introns_on_the_minus_strand():
    found = entries()
    exons = ((101, 120), (126, 155), (161, 170))
    assert found["YAA002C_id001"].sequence == piece("chrI", *exons, minus=True)
    five = found["YAA002C_id001:five_prime_UTR"]
    assert five.sequence == piece("chrI", (151, 155), (161, 170), minus=True)
    assert [(e.start, e.stop) for e in five.regions[0].exons] == [
        (151, 155),
        (161, 170),
    ]
    assert found["YAA002C_id001:three_prime_UTR"].sequence == piece(
        "chrI", (101, 110), minus=True
    )


def test_links_utrs_and_mrnas_both_ways():
    found = entries()
    mrna = found["YAA001W_id001"]
    assert {(r.sequence_id, r.relationship) for r in mrna.related_sequences} == {
        ("SGD_CODING:YAA001W_id001:five_prime_UTR", "five_prime_utr"),
        ("SGD_CODING:YAA001W_id001:three_prime_UTR", "three_prime_utr"),
    }
    utr = found["YAA001W_id001:three_prime_UTR"]
    assert [(r.sequence_id, r.relationship) for r in utr.related_sequences] == [
        ("SGD_CODING:YAA001W_id001", "mrna")
    ]


def test_skips_dubious_orfs_and_short_utrs():
    found = entries()
    assert not any(a.startswith("YAA003W") for a in found)
    assert "Q0001_mRNA" in found
    assert "Q0001_mRNA:five_prime_UTR" not in found


def test_uses_rnacentral_assembly_and_chromosome_names():
    found = entries()
    region = found["Q0001_mRNA"].regions[0]
    assert (region.assembly_id, region.chromosome) == ("R64-1-1", "Mito")
    mrna = found["YAA001W_id001"]
    assert mrna.regions[0].chromosome == "I"
    assert mrna.database == "SGD_CODING"
    assert mrna.accession == "SGD_CODING:YAA001W_id001"
    assert mrna.ncbi_tax_id == 559292
    assert mrna.gene == "AAA1"
    assert mrna.locus_tag == "YAA001W"
    assert mrna.url == "https://www.yeastgenome.org/locus/S000000001"


def test_handles_an_intron_at_the_end_of_the_mrna():
    found = entries()
    mrna = found["YAA005C_mRNA"]
    assert mrna.sequence == piece("chrmt", (45, 55), minus=True)
    assert [(e.start, e.stop) for e in mrna.regions[0].exons] == [(45, 55)]
