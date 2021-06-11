# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import attr
import pytest
from io import open

from rnacentral_pipeline.databases.data import Exon
from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.data import SequenceRegion
from rnacentral_pipeline.databases.data import CoordinateSystem

from rnacentral_pipeline.databases.rgd import helpers as rgd


@pytest.fixture(scope="module")  # pylint: disable=no-member
def simple_entry():
    with open("data/rgd/rat_ncrna.tsv", "r") as raw:
        return next(rgd.as_rows(raw))


@pytest.fixture(scope="module")  # pylint: disable=no-member
def tricky_entry():
    with open("data/rgd/rat_ncrna.tsv", "r") as raw:
        for entry in rgd.as_rows(raw):
            if entry["GENE_RGD_ID"] == "10401674":
                return entry
    raise ValueError("What!")


@pytest.fixture(scope="module")  # pylint: disable=no-member
def rat_ncrna():
    with open("data/rgd/rat_ncrna.tsv", "r") as raw:
        return list(rgd.as_rows(raw))


@pytest.fixture(scope="module")  # pylint: disable=no-member
def rat_protein():
    with open("data/rgd/rat_protein.tsv", "r") as raw:
        return list(rgd.as_rows(raw))


@pytest.fixture(scope="module")  # pylint: disable=no-member
def rat_multi_locus():
    with open("data/rgd/multi-locus.tsv", "r") as raw:
        return list(rgd.as_rows(raw))


@pytest.fixture(scope="module")  # pylint: disable=no-member
def sequences():
    with rgd.indexed("data/rgd/sequences.fa.gz") as indexed:
        yield indexed


def test_can_fetch_all_organisms():
    assert rgd.known_organisms() == ["rat"]


def test_can_fetch_accession_from_entry(simple_entry):
    assert rgd.accession(simple_entry) == "RRID:RGD_5687330"


def test_can_fetch_accession_from_entry_with_index(simple_entry):
    assert rgd.accession(simple_entry, 0) == "RRID:RGD_5687330:0"


def test_can_determine_taxid(simple_entry):
    assert rgd.taxid(simple_entry) == 10116


# @pytest.mark.parametrize('entry', rat_ncrna)
def test_can_detect_if_is_ncrna(rat_ncrna):
    assert rgd.is_ncrna(rat_ncrna[0]) is True


# @pytest.mark.parametrize('entry', ())
def test_can_detect_if_not_ncrna(rat_protein):
    assert rgd.is_ncrna(rat_protein[0]) is False


def test_can_generate_url(simple_entry):
    assert (
        rgd.url(simple_entry)
        == "https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=5687330"
    )


def test_can_generate_exons(simple_entry):
    assert rgd.regions(simple_entry) == [
        SequenceRegion(
            chromosome="14",
            strand=1,
            exons=[Exon(start=46643222, stop=46645093)],
            assembly_id="",
            coordinate_system=CoordinateSystem.zero_based(),
        )
    ]


def test_can_generate_reasonable_description(simple_entry):
    assert (
        rgd.description(simple_entry) == "Rattus norvegicus 18S ribosomal RNA (Rn18s)"
    )


def test_can_build_xrefs(simple_entry):
    assert rgd.xref_data(simple_entry) == {
        "genbank": ["AABR07015078", "AH001747", "DQ623540", "NR_046237", "V01270"],
        "ncbi_gene": ["100861533"],
    }


@pytest.mark.xfail()
def test_can_fetch_sequence(simple_entry):
    with rgd.indexed("data/rgd/sequences.fa.gz") as indexed:
        data = rgd.sequences_for(simple_entry, indexed)
        seqs = [d[0] for d in data]
        # exons = [d[1] for d in data]
        assert seqs == [
            "TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGTGGGGGGGACGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCCTCCCGGCTCCGGCCGGGGGTCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCTCCGTGGCGGCGACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGGGTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGTCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACTGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA"
        ]
        # assert exons == []


@pytest.mark.skip()
def test_fails_without_existing_sequence(simple_entry, sequences):
    entry = dict(simple_entry)
    entry["GENE_RGD_ID"] = "something-made-up"
    with pytest.raises(Exception):
        rgd.sequences_for(entry, sequences)


@pytest.mark.xfail()
def test_can_build_correct_entry(simple_entry, sequences):
    entries = rgd.as_entries(simple_entry, sequences)
    assert len(entries) == 1
    assert attr.asdict(entries[0]) == attr.asdict(
        Entry(
            primary_id="5687330",
            accession="RRID:RGD_5687330",
            ncbi_tax_id=10116,
            database="RGD",
            sequence="TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGTGGGGGGGACGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCCTCCCGGCTCCGGCCGGGGGTCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCTCCGTGGCGGCGACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGGGTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGTCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACTGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA",
            exons=[],
            rna_type="rRNA",
            url="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=5687330",
            seq_version="1",
            xref_data={
                "genbank": [
                    "AABR07015078",
                    "AH001747",
                    "DQ623540",
                    "NR_046237",
                    "V01270",
                ],
                "ncbi_gene": ["100861533"],
            },
            species="Rattus norvegicus",
            common_name="Norway rat",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; "
                "Muridae; Murinae; Rattus; Rattus norvegicus"
            ),
            gene="Rn18s",
            locus_tag="Rn18s",
            description="Rattus norvegicus 18S ribosomal RNA (Rn18s)",
            references=[
                Reference(
                    authors="Shimoyama M, De Pons J, Hayman GT, Laulederkind SJ, Liu W, Nigam R, Petri V, Smith JR, Tutaj M, Wang SJ, Worthey E, Dwinell M, Jacob H.",
                    location="Nucleic Acids Res 43(database issue):D743-50 (2015)",
                    title="The Rat Genome Database 2015: genomic, phenotypic and environmental variations and disease",
                    pmid=25355511,
                    doi=u"10.1093/nar/gku1026",
                ),
                Reference(
                    authors="Girard A, Sachidanandam R, Hannon GJ, Carmell MA.",
                    location="Nature 442(7099):199-202 (2006)",
                    title="A germline-specific class of small RNAs binds mammalian Piwi proteins",
                    pmid=16751776,
                    doi=u"10.1038/nature04917",
                ),
                Reference(
                    authors="Choi YC.",
                    location="J Biol Chem 260(23):12769-12772 (1985)",
                    title="Structural organization of ribosomal RNAs from Novikoff hepatoma. I. Characterization of fragmentation products from 40 S subunit",
                    pmid=3930503,
                    doi=None,
                ),
                Reference(
                    authors="Subrahmanyam CS, Cassidy B, Busch H, Rothblum LI.",
                    location="Nucleic Acids Res 10(12):3667-3680 (1982)",
                    title="Nucleotide sequence of the region between the 18S rRNA sequence and the 28S rRNA sequence of rat ribosomal DNA",
                    pmid=6287418,
                    doi=u"10.1093/nar/10.12.3667",
                ),
                Reference(
                    authors="Rothblum LI, Reddy R, Cassidy B.",
                    location="Nucleic Acids Res 10(22):7345-7362 (1982)",
                    title="Transcription initiation site of rat ribosomal DNA",
                    pmid=6296773,
                    doi=u"10.1093/nar/10.22.7345",
                ),
                Reference(
                    authors="Chan YL, Olvera J, Wool IG.",
                    location="Nucleic Acids Res 11(22):7819-7831 (1983)",
                    title="The structure of rat 28S ribosomal ribonucleic acid inferred from the sequence of nucleotides in a gene",
                    pmid=6316273,
                    doi=u"10.1093/nar/11.22.7819",
                ),
                Reference(
                    authors="Chan YL, Gutell R, Noller HF, Wool IG.",
                    location="J Biol Chem 259(1):224-230 (1984)",
                    title="The nucleotide sequence of a rat 18 S ribosomal ribonucleic acid gene and a proposal for the secondary structure of 18 S ribosomal ribonucleic acid",
                    pmid=6323401,
                    doi=None,
                ),
                Reference(
                    authors="Hadjiolov AA, Georgiev OI, Nosikov VV, Yavachev LP.",
                    location="Nucleic Acids Res 12(8):3677-3693 (1984)",
                    title="Primary and secondary structure of rat 28 S ribosomal RNA",
                    pmid=6328433,
                    doi=u"10.1093/nar/12.8.3677",
                ),
            ],
        )
    )


@pytest.mark.xfail()
def test_produces_none_for_no_sequence(tricky_entry, sequences):
    assert attr.asdict(rgd.as_entries(tricky_entry, sequences)) == attr.asdict(
        Entry(
            primary_id="10401674",
            accession="RRID:RGD_10401674",
            ncbi_tax_id=10116,
            database="RGD",
            sequence="",
            exons=[],
            rna_type="ncrna",
            url="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=10401674",
            seq_version="1",
            xref_data={},
            gene="Carmn",
            locus_tag="Carmn",
            gene_synonyms=["Mir143hg"],
            description="Rattus norvegicus cardiac mesoderm enhancer-associated non-coding RNA (Carmn)",
        )
    )


@pytest.mark.xfail()
def test_can_create_entries_for_all_ncrna(rat_ncrna, sequences):
    entries = [rgd.as_entries(e, sequences) for e in rat_ncrna]
    assert len(entries) == 15
    assert len([e for e in entries if e]) == 14


@pytest.mark.xfail()
def test_can_handle_having_multiple_locations(rat_multi_locus, sequences):
    entries = []
    for raw in rat_multi_locus:
        entries.extend(rgd.as_entries(raw, sequences))

    assert len(entries) == 3
    assert attr.asdict(entries[0]) == attr.asdict(
        Entry(
            primary_id="7706003",
            accession="RRID:RGD_7706003:1",
            ncbi_tax_id=10116,
            database="RGD",
            sequence="TCTATGCTGTGATCTTAAAGCAGACAGCCACTCAAGGCCCTAGCAGGATCGGCACTTGTGTCTCCTTTCGGGGGTGCCAGCCTTCTCGCCTTGAGTTTGTTTTCTTCCTGAAACTCTCCAGTCATCACAAGCCTATGTGACCAGATGCTCAGCCCCTCTGAGTTCTGGCATCTACTACTTCTGTTTCTTCAGCTCCATCCAATCCTGTTGGTCCCAGTCCTGGCCTTGACTTTGGCTACCAAGCACTCACCTGGCGTCTACCTCTCCTGAACCCGGCTTCTTAAGGATACCGAAGTATATGGAGGGGGGCACAGCTGTCTGCAACTCCAGCTCAGGTGATCTGACACCCTCTTTTGGCCTCCACGGACAGCTAAACGCACACAGTTGCAGAACGACCCAGATCATCGGTGGATATTCAGCATGACCTCATCAATGCAAGCAGTGCCCAGACGAAGAATGAGGACATCAGCAGCCCCACCCCCACCCCCTTGGTATCTCCCAGCCTCTTTCAACCACTGGCCAAACAAGGAGCCTCTATCCTGACCTTTAACATCGAAGATTGGTTTGGGCTGGTTTTGAGTTGACATAAATGGGATTGCCCAGTACACATTCCTGTCTGGATTCTTTGCTCCGCCTTATGTTTTCCAGCTTGGCCTGTGTTTTGCAGAGCAAACCATTCATTTTCCACAGTGCAGTAGAGCGTTCTACACTTCCAGAGTGTGGAACACAATTAAGAGTTCATGGAGAAAGAGAAGAGCTGGAGAGTGAGCGAG",
            exons=[],
            rna_type="other",
            url="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=7706003",
            seq_version="1",
            xref_data={
                "ncbi_gene": ["102549948"],
                "ensembl": ["ENSRNOG00000052150"],
                "genbank": [
                    "AABR07067034",
                    "AAHX01058687",
                    "AAHX01058688",
                    "XR_001839772",
                    "XR_348953",
                    "XR_356918",
                    "XR_356920",
                ],
            },
            species="Rattus norvegicus",
            common_name="Norway rat",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; "
                "Muridae; Murinae; Rattus; Rattus norvegicus"
            ),
            gene="LOC102549948",
            locus_tag="LOC102549948",
            description="Rattus norvegicus uncharacterized LOC102549948",
            references=[
                Reference(
                    authors="Shimoyama M, De Pons J, Hayman GT, Laulederkind SJ, Liu W, Nigam R, Petri V, Smith JR, Tutaj M, Wang SJ, Worthey E, Dwinell M, Jacob H.",
                    location="Nucleic Acids Res 43(database issue):D743-50 (2015)",
                    title="The Rat Genome Database 2015: genomic, phenotypic and environmental variations and disease",
                    pmid=25355511,
                    doi=u"10.1093/nar/gku1026",
                )
            ],
        )
    )

    assert attr.asdict(entries[1]) == attr.asdict(
        Entry(
            primary_id="7706003",
            accession="RRID:RGD_7706003:2",
            ncbi_tax_id=10116,
            database="RGD",
            sequence="TTGACTTGGATCATCTACTCATCGATCTGTAAACGTGCATCAAGATGTCAGATAACTACTCATGGAAATTGGTTATGCTCCCAGACAAAGAAACCAAAATAAACCAGGCCTTTGATCAGCAAGCAGTCAGCACGACCCAGAGGCAGAGAGATACCGAAGTATATGGAGGGGGGCACAGCTGTCTGCAACTCCAGCTCAGGTGATCTGACACCCTCTTTTGGCCTCCACGGACAGCTAAACGCACACAGTTGCAGAACGACCCAGATCATCGGTGGATATTCAGCATGACCTCATCAATGCAAGCAGTGCCCAGACGAAGAATGAGGACATCAGCAGCCCCACCCCCACCCCCTTGGTATCTCCCAGCCTCTTTCAACCACTGGCCAAACAAGGAGCCTCTATCCTGACCTTTAACATCGAAGATTGGTTTGGGCTGGTTTTGAGTTGACATAAATGGGATTGCCCAGTACACATTCCTGTCTGGATTCTTTGCTCCGCCTTATGTTTTCCAGCTTGGCCTGTGTTTTGCAGAGCAAACCATTCATTTTCCACAGTGCAGTAGAGCGTTCTACACTTCCAGAGTGTGGAACACAATTAAGAGTTCATGGAGAAAGAGAAGAGCTGGAGAGTGAGCGAG",
            exons=[],
            rna_type="other",
            url="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=7706003",
            seq_version="1",
            xref_data={
                "ncbi_gene": ["102549948"],
                "ensembl": ["ENSRNOG00000052150"],
                "genbank": [
                    "AABR07067034",
                    "AAHX01058687",
                    "AAHX01058688",
                    "XR_001839772",
                    "XR_348953",
                    "XR_356918",
                    "XR_356920",
                ],
            },
            species="Rattus norvegicus",
            common_name="Norway rat",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; "
                "Muridae; Murinae; Rattus; Rattus norvegicus"
            ),
            gene="LOC102549948",
            locus_tag="LOC102549948",
            description="Rattus norvegicus uncharacterized LOC102549948",
            references=[
                Reference(
                    authors="Shimoyama M, De Pons J, Hayman GT, Laulederkind SJ, Liu W, Nigam R, Petri V, Smith JR, Tutaj M, Wang SJ, Worthey E, Dwinell M, Jacob H.",
                    location="Nucleic Acids Res 43(database issue):D743-50 (2015)",
                    title="The Rat Genome Database 2015: genomic, phenotypic and environmental variations and disease",
                    pmid=25355511,
                    doi=u"10.1093/nar/gku1026",
                )
            ],
        )
    )

    assert attr.asdict(entries[2]) == attr.asdict(
        Entry(
            primary_id="7706003",
            accession="RRID:RGD_7706003:3",
            ncbi_tax_id=10116,
            database="RGD",
            sequence="ACTAAGCTACTAACACAGAGAGATAATTCATTCTTAGATGGGATGTAATGTACTGCATCTGTTAGGAGGCACTTTCCGAATTAATTACTCCGTCACTCCCTCCCCTGCCCCATCTCTTTCCTTAGTTCTTACAGTCCTAGTCGTGGCTCACTAGCCTGACCCTGCCTGAAGAGATAGGAAAGGAACTGTGCTCCACGTAAATGAACCTGCCTCACCCTCAAAATAATCCAGAGCCTGGATCATCCCACATGGACGCTACTGGACCGGTGGCCTGCAAGATAGAACCAAGAAACTTTCTTCCCAACTGTTATAAAAAAGAACAGAACAGAAGAGAACAAAACACAGCAATCCTTAGATTCTCTCCCAGCCCACTTCCCCTCTGGAAACCCATTTACTGTCTCCTTCTAGTGTATTTCAATAAATTCTCCTGTGTTCCAGTTAA",
            exons=[],
            rna_type="other",
            url="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=7706003",
            seq_version="1",
            xref_data={
                "ncbi_gene": ["102549948"],
                "ensembl": ["ENSRNOG00000052150"],
                "genbank": [
                    "AABR07067034",
                    "AAHX01058687",
                    "AAHX01058688",
                    "XR_001839772",
                    "XR_348953",
                    "XR_356918",
                    "XR_356920",
                ],
            },
            species="Rattus norvegicus",
            common_name="Norway rat",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; "
                "Muridae; Murinae; Rattus; Rattus norvegicus"
            ),
            gene="LOC102549948",
            locus_tag="LOC102549948",
            description="Rattus norvegicus uncharacterized LOC102549948",
            references=[
                Reference(
                    authors="Shimoyama M, De Pons J, Hayman GT, Laulederkind SJ, Liu W, Nigam R, Petri V, Smith JR, Tutaj M, Wang SJ, Worthey E, Dwinell M, Jacob H.",
                    location="Nucleic Acids Res 43(database issue):D743-50 (2015)",
                    title="The Rat Genome Database 2015: genomic, phenotypic and environmental variations and disease",
                    pmid=25355511,
                    doi=u"10.1093/nar/gku1026",
                )
            ],
        )
    )
