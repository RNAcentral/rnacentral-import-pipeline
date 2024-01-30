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

import os
from pathlib import Path

import attr
import pytest

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.ena import context
from rnacentral_pipeline.databases.ena import mapping as tpa
from rnacentral_pipeline.databases.ena import parser
from rnacentral_pipeline.databases.helpers.hashes import md5


def parse(path):
    builder = context.ContextBuilder()
    builder.with_dr(path)
    ctx = builder.context()
    return parser.parse(ctx, path)


@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/ena/tpa/lncrnadb/mapping.tsv", 62),
        ("data/ena/tpa/snopy/mapping.tsv", 2634),
        ("data/ena/tpa/srpdb/mapping.tsv", 855),
        ("data/ena/tpa/tmrna/mapping.tsv", 21318),
        ("data/ena/tpa/wormbase/mapping.tsv", 27665),
        ("data/ena/tpa/tair/mapping.tsv", 1290),
    ],
)
def test_can_parse_complete_tpa_files(filename, count):
    with open(filename, "r") as raw:
        loaded = list(tpa.parse_tpa_file(raw))
        assert len(loaded) == count


@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/ena/tpa/lncrnadb/mapping.tsv", 62),
        ("data/ena/tpa/snopy/mapping.tsv", 2634),
        ("data/ena/tpa/srpdb/mapping.tsv", 855),
        ("data/ena/tpa/tmrna/mapping.tsv", 21318),
        ("data/ena/tpa/wormbase/mapping.tsv", 27665),
        # ('data/ena/tpa/tair/mapping.tsv', 1290),
    ],
)
def test_can_produce_correct_number_of_tpa_keys_from_tpa_file(filename, count):
    with open(filename, "r") as raw:
        loaded = list(tpa.parse_tpa_file(raw))
        keys = set(tpa.tpa_key(t) for t in loaded)
        assert len(keys) == count


def test_can_build_correct_lncrnadb_tpas():
    with open("data/ena/tpa/lncrnadb/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert data == tpa.GenericTpa(
        "LNCRNADB",
        "101",
        "Kcnq1ot1",
        "HG975405",
        None,
    )


def test_can_build_correct_snopy_tpas():
    with open("data/ena/tpa/snopy/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert data == tpa.GenericTpa(
        "SNOPY",
        "Arabidopsis_thaliana300001",
        "SnoR1b",
        "LN809305",
        None,
    )


def test_can_build_correct_srpdb_tpas():
    with open("data/ena/tpa/srpdb/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert data == tpa.GenericTpa(
        "SRPDB",
        "Acin.baum._CP000521",
        None,
        "HG323367",
        None,
    )


def test_can_build_correct_tmrna_tpas():
    with open("data/ena/tpa/tmrna/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert data == tpa.GenericTpa(
        "TMRNA_WEB",
        "Acary_marin_MBIC11",
        None,
        "HG525190",
        None,
    )


def test_can_build_correct_wormbase_tpas():
    with open("data/ena/tpa/wormbase/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert data == tpa.GenericTpa(
        "WORMBASE",
        "WBGene00001734",
        "ZK643.8b",
        "BX284603",
        None,
    )


def test_can_build_correct_tair_tpas():
    with open("data/ena/tpa/tair/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert data == tpa.GenericTpa(
        "TAIR",
        "1000429427",
        "AT5G52495",
        "CP002688",
        None,
    )


def transform_first(db_name):
    embl_file = os.path.join("data", "ena", "tpa", db_name, "entry.embl")
    embl = Path(embl_file)
    entry = next(parse(embl))

    tpa_file = os.path.join("data", "ena", "tpa", db_name, "mapping.tsv")
    with open(tpa_file, "r") as raw:
        tpa_entry = next(tpa.parse_tpa_file(raw))

    urls = tpa.UrlBuilder()
    updated = tpa_entry.transform(entry)
    return urls.transform(updated)


def test_can_transform_correct_lncrnadb_entry():
    transformed = attr.asdict(transform_first("lncrnadb"))
    assert (
        md5(transformed["sequence"].encode("utf-8"))
        == "3b9990f42f263eb0894f92fada30273d"
    )
    assert len(transformed["sequence"]) == 32753

    result = attr.asdict(
        Entry(
            primary_id="101",
            accession="HG975405.1:1..32753:ncRNA:LNCRNADB:101",
            ncbi_tax_id=9606,
            database="LNCRNADB",
            sequence="",
            regions=[],
            rna_type="lncRNA",
            url="http://www.lncrnadb.org/Detail.aspx?TKeyID=101",
            seq_version="1",
            note_data={
                "ontology": ["ECO:0000305", "SO:0001903"],
                "text": [
                    "biotype:sense_intronic",
                    "transcribed from intron 10 of the maternally expressed Kncq1 (KvLQT1) gene",
                ],
            },
            xref_data={
                "lncrnadb": ["101", "Kcnq1ot1"],
                "ena_refs": {
                    "LNCRNADB": ("101", "Kcnq1ot1"),
                },
            },
            optional_id="Kcnq1ot1",
            non_coding_id="HG975405.1:1..32753:ncRNA",
            is_composite="Y",
            species="Homo sapiens",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;"
                " Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; "
                "Catarrhini; Hominidae; Homo; Homo sapiens"
            ),
            common_name="human",
            project="PRJEB6238",
            product="Long non-coding sense-intronic RNA Kcnq1ot1",
            description="Homo sapiens (human) Long non-coding sense-intronic RNA Kcnq1ot1",
            keywords="RNAcentral; TPA; TPA:specialist_db",
            parent_accession="HG975405",
            gene="Kncq1",
            mol_type="transcribed RNA",
            experiment=(
                "EXISTENCE: lncRNAdb literature review [PMID: "
                "15340049,18299392,18951091,18848501,19144718,17917697, "
                "20573698,10393948,15516932,17242189,10369866,11813134, "
                "16702402,16575194,16965397,10958646,21345374,15590939, "
                "15459184,21172659,21576366,22406755]"
            ),
        )
    )

    # We remove sequence because it is too big to really look at
    del transformed["sequence"]
    del result["sequence"]

    # References are too long to compare either so I check the count
    assert len(transformed["references"]) == 24
    del transformed["references"]
    del result["references"]
    assert transformed == result


def test_can_transform_correct_srpdb_entry():
    transformed = attr.asdict(transform_first("srpdb"))
    result = attr.asdict(
        Entry(
            primary_id="Acin.baum._CP000521",
            accession="HG323367.1:1..116:ncRNA:SRPDB:Acin.baum._CP000521",
            ncbi_tax_id=400667,
            database="SRPDB",
            sequence="GGTAGCCCCGCTGGTGGCTTCGCAATTGTACTCTGTGAACCCCGCCAGGACCGGAAGGTAGCAACGGTAGCAGATCTATGATGTGCCGAAGTTTTGCTAGTGGGGTTGCCACCATT",
            regions=[],
            rna_type="SRP_RNA",
            url="http://rnp.uthscsa.edu/rnp/SRPDB/rna/sequences/fasta/Acin.baum._CP000521",
            seq_version="1",
            xref_data={
                "ena_refs": {
                    "SRPDB": ("Acin.baum._CP000521", None),
                }
            },
            note_data={
                "ontology": ["ECO:0000305", "GO:0006617", "GO:0048501", "SO:0000590"],
                "text": ["alignment group:Small 4.5S Bacteria (SB)"],
            },
            species="Acinetobacter baumannii ATCC 17978",
            lineage=(
                "Bacteria; Proteobacteria; Gammaproteobacteria; "
                "Pseudomonadales; Moraxellaceae; Acinetobacter; "
                "Acinetobacter calcoaceticus/baumannii complex; "
                "Acinetobacter baumannii ATCC 17978"
            ),
            keywords="RNAcentral; TPA; TPA:specialist_db",
            description="Acinetobacter baumannii ATCC 17978 signal recognition particle RNA",
            mol_type="transcribed RNA",
            gene="SRP RNA",
            non_coding_id="HG323367.1:1..116:ncRNA",
            is_composite="Y",
            project="PRJEB4384",
            parent_accession="HG323367",
            product="signal recognition particle RNA",
        )
    )
    assert len(transformed["references"]) == 3
    del transformed["references"]
    del result["references"]

    assert transformed == result


def test_can_transform_correct_snopy_entry():
    transformed = attr.asdict(transform_first("snopy"))
    result = attr.asdict(
        Entry(
            primary_id="Arabidopsis_thaliana300001",
            accession="LN809305.1:1..93:ncRNA:SNOPY:Arabidopsis_thaliana300001",
            ncbi_tax_id=3702,
            database="SNOPY",
            sequence="GGCGAGGATGAATAATGCTAAATTTCTGACACCTCTTGTATGAGGAGAGATTGATAACCTCTCCTTTGAGCACATTATGCAATACTCTGAGCC",
            regions=[],
            rna_type="snoRNA",
            url="http://snoopy.med.miyazaki-u.ac.jp/snorna_db.cgi?mode=sno_info&id=Arabidopsis_thaliana300001",
            seq_version="1",
            note_data={
                "ontology": ["ECO:0000305", "GO:0005730", "GO:0006396", "SO:0000275"],
                "text": [
                    "Modification:C/D",
                    "Target:25S rRNA, 18S rRNA",
                    "Organisation:Poly",
                ],
            },
            description="Arabidopsis thaliana (thale cress) small nucleolar RNA SnoR1b",
            species="Arabidopsis thaliana",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; "
                "Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; "
                "Gunneridae; Pentapetalae; rosids; malvids; Brassicales; "
                "Brassicaceae; Camelineae; Arabidopsis; Arabidopsis thaliana"
            ),
            common_name="thale cress",
            keywords="RNAcentral; TPA; TPA:specialist_db",
            mol_type="genomic DNA",
            gene="SnoR1b",
            is_composite="Y",
            non_coding_id="LN809305.1:1..93:ncRNA",
            # experiment='EXISTENCE:Curator Inference ECO0000305',
            product="small nucleolar RNA SnoR1b",
            project="PRJEB8122",
            parent_accession="LN809305",
        )
    )

    assert len(transformed["references"]) == 2
    del transformed["references"]
    del result["references"]

    assert transformed == result


def test_can_transform_correct_wormbase_entry():
    transformed = attr.asdict(transform_first("wormbase"))
    result = attr.asdict(
        Entry(
            primary_id="WBGene00001734",
            accession="BX284603.4:8962295..8965569:misc_RNA:WORMBASE:WBGene00001734",
            ncbi_tax_id=6239,
            database="WORMBASE",
            sequence="",
            regions=[],
            rna_type="misc_RNA",
            url="http://www.wormbase.org/species/c_elegans/gene/WBGene00001734",
            seq_version="4",
            xref_data={
                "ena_refs": {
                    "BIOSAMPLE": ("SAMEA3138177", None),
                    "WORMBASE": ("WBGene00001734", "ZK643.8b"),
                },
            },
            description="Caenorhabditis elegans Non-coding transcript of protein-coding gene grl-25",
            species="Caenorhabditis elegans",
            lineage=(
                "Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida;"
                " Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis; "
                "Caenorhabditis elegans"
            ),
            mol_type="genomic DNA",
            chromosome="III",
            locus_tag="CELE_ZK643.8",
            gene="grl-25",
            optional_id="ZK643.8b",
            is_composite="Y",
            non_coding_id="BX284603.4:8962295..8965569:misc_RNA",
            product="Non-coding transcript of protein-coding gene grl-25",
            project="PRJNA13758",
            parent_accession="BX284603",
            standard_name="ZK643.8b",
        )
    )

    assert len(transformed["references"]) == 3
    assert len(transformed["regions"]) == 0
    assert len(transformed["sequence"]) == 2767
    del transformed["references"]
    del result["references"]
    del transformed["sequence"]
    del result["sequence"]
    del transformed["regions"]
    del result["regions"]
    assert transformed == result


def test_build_correct_tpa_key_for_snopy_entries():
    raw = Path("data/ena/tpa/snopy/entry.embl")
    entry = next(parse(raw))
    assert tpa.tpa_key(entry) == ("LN809305", None)


def test_builds_correct_tpa_key_for_snopy_tpa():
    with open("data/ena/tpa/snopy/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert tpa.tpa_key(data) == ("LN809305", None)


def test_builds_correct_tpa_key_for_wormbase_tpa():
    with open("data/ena/tpa/wormbase/mapping.tsv", "r") as raw:
        data = next(tpa.parse_tpa_file(raw))
    assert tpa.tpa_key(data) == ("BX284603", "ZK643.8b")


def test_knows_if_it_has_mappings():
    mapping = tpa.load_file("data/ena/tpa/snopy/mapping.tsv")
    raw = Path("data/ena/tpa/snopy/entry.embl")
    entry = next(parse(raw))
    assert mapping.has_tpa_for(entry) is True


def test_knows_if_does_not_have_entry_for():
    mapping = tpa.load_file("data/ena/tpa/snopy/mapping.tsv")
    raw = Path("data/ena/tpa/lncrnadb/entry.embl")
    entry = next(parse(raw))
    assert mapping.has_tpa_for(entry) is False


def test_can_fetch_tpa_for_entry():
    mapping = tpa.load_file("data/ena/tpa/snopy/mapping.tsv")
    raw = Path("data/ena/tpa/snopy/entry.embl")
    entry = next(parse(raw))
    tpas = list(mapping.find_tpas(entry))
    assert len(tpas) == 1
    assert tpas[0] == tpa.GenericTpa(
        "SNOPY",
        "Arabidopsis_thaliana300001",
        "SnoR1b",
        "LN809305",
        None,
    )


def test_fetch_tpa_for_non_exist_returns_empty_list():
    mapping = tpa.load_file("data/ena/tpa/snopy/mapping.tsv")
    raw = Path("data/ena/tpa/lncrnadb/entry.embl")
    entry = next(parse(raw))
    assert list(mapping.find_tpas(entry)) == []


def test_can_map_snopy_entries():
    mapping = tpa.load_file("data/ena/tpa/snopy/mapping.tsv")

    raw = Path("data/ena/tpa/snopy/entry.embl")
    entries = list(parse(raw))

    mapped = list(tpa.apply(mapping, entries))
    assert len(mapped) == 1
    assert mapped[0].database == "SNOPY"
    assert (
        mapped[0].accession == "LN809305.1:1..93:ncRNA:SNOPY:Arabidopsis_thaliana300001"
    )


def test_can_will_not_alter_entries_from_other_dbs():
    mapping = tpa.load_file("data/ena/tpa/snopy/mapping.tsv")

    raw = Path("data/ena/tpa/lncrnadb/entry.embl")
    entries = list(parse(raw))

    mapped = list(tpa.apply(mapping, entries))
    assert len(mapped) == 1
    assert mapped[0].database == "ENA"
    assert mapped[0].accession == "HG975405.1:1..32753:ncRNA"


def test_can_apply_wormbase_tpas():
    mapping = tpa.load_file("data/ena/tpa/wormbase/mapping.tsv")

    raw = Path("data/ena/tpa/wormbase/entry.embl")
    entries = list(parse(raw))
    assert entries

    mapped = list(tpa.apply(mapping, entries))
    assert len(mapped) == len(entries)
    assert mapped[0].database == "WORMBASE"
    assert (
        mapped[0].accession
        == "BX284603.4:8962295..8965569:misc_RNA:WORMBASE:WBGene00001734"
    )


def test_can_apply_mirbase_tpas():
    mapping = tpa.load_file("data/ena/tpa/mirbase/mapping.tsv")
    raw = Path("data/ena/tpa/mirbase/entry.embl")
    entries = list(parse(raw))
    assert entries
    mapped = list(tpa.apply(mapping, entries))
    assert len(mapped) == len(entries)

    assert mapped[0].database == "MIRBASE"
    assert mapped[0].accession == "LM611181.1:1..180:precursor_RNA:MIRBASE:MI0016048"
    assert mapped[0].optional_id == "hsa-mir-3648-1"


@pytest.mark.parametrize("name", [tpa.internal_database_name(n) for n in tpa.DATABASES])
def test_url_builder_can_build_url(name):
    assert hasattr(tpa.UrlBuilder(), name.lower())


@pytest.mark.parametrize(
    "filename,fails",
    [
        ("data/ena/tpa/lncrnadb/mapping.tsv", True),
        ("data/ena/tpa/mirbase/mapping.tsv", True),
        ("data/ena/tpa/snopy/mapping.tsv", True),
        ("data/ena/tpa/srpdb/mapping.tsv", True),
        ("data/ena/tpa/tair/mapping.tsv", True),
        ("data/ena/tpa/tmrna/mapping.tsv", True),
        ("data/ena/tpa/wormbase/mapping.tsv", True),
        pytest.param(
            "data/ena/tpa/combined/mapping.tsv",
            False,
            marks=pytest.mark.skip(reason="SGD is broken"),
        ),
    ],
)
def test_fails_to_parse_if_no_tpas(filename, fails):
    if fails:
        with pytest.raises(Exception):
            tpa.load_file(filename).validate()
    else:
        assert tpa.load_file(filename).validate() is True
