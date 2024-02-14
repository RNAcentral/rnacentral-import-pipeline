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

from rnacentral_pipeline.databases import data as dat

from .helpers import entries_for, entry_for, has_entry_for, parse_with_family


@pytest.fixture(scope="module")  # pylint: disable=no-member
def human_1():
    return parse_with_family(
        "data/ensembl/Homo_sapiens.GRCh38.chromosome.1.dat",
        gencode_file="data/gencode/human-transcripts.gff3",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def human_12():
    with open("data/ensembl/excluded.txt", "r") as ex:
        return parse_with_family(
            "data/ensembl/Homo_sapiens.GRCh38.chromosome.12.dat",
            gencode_file="data/gencode/human-transcripts.gff3",
            excluded_file=ex,
        )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def human_x():
    return parse_with_family("data/ensembl/Homo_sapiens.GRCh38.chromosome.X.dat")


@pytest.fixture(scope="module")  # pylint: disable=no-member
def macaca():
    return parse_with_family(
        "data/ensembl/Macaca_mulatta.Mmul_8.0.1.chromosome.1.dat",
        gencode_file="data/gencode/human-transcripts.gff3",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def mouse_3():
    return parse_with_family("data/ensembl/Mus_musculus.GRCm38.chromosome.3.dat")


@pytest.fixture(scope="module")  # pylint: disable=no-member
def cow_8():
    return parse_with_family(
        "data/ensembl/Bos_taurus.ARS-UCD1.2.primary_assembly.8.dat"
    )


@pytest.mark.slow
def test_it_sets_primary_id_to_versionless_transcript_id(human_12):
    assert entry_for(human_12, "ENST00000516089.1").primary_id == "ENST00000516089"


@pytest.mark.slow
def test_it_generates_correct_seq_version(human_12):
    assert entry_for(human_12, "ENST00000516089.1").seq_version == "1"


@pytest.mark.slow
def test_sets_optional_id_to_gene_id(human_12):
    assert entry_for(human_12, "ENST00000516089.1").optional_id == "ENSG00000251898.1"


@pytest.mark.slow
def test_it_gets_gene_id_to_locus(human_12):
    assert entry_for(human_12, "ENST00000516089.1").gene == "SCARNA11"


@pytest.mark.slow
def test_it_gets_the_locus_tag(human_12):
    assert entry_for(human_12, "ENST00000516089.1").locus_tag == "SCARNA11"


@pytest.mark.slow
def test_it_sets_rna_type_to_snRNA(human_12):
    assert entry_for(human_12, "ENST00000516089.1").rna_type == "SO:0002095"
    assert entry_for(human_12, "ENST00000540226.1").rna_type == "SO:0001877"


@pytest.mark.slow
def test_it_sets_product_to_scaRNA(human_12):
    assert entry_for(human_12, "ENST00000516089.1").product == "scaRNA"
    assert entry_for(human_12, "ENST00000516089.1").rna_type == "SO:0002095"


@pytest.mark.slow
def test_it_sets_accession_to_transcript_id(human_12):
    assert entry_for(human_12, "ENST00000540868.1").accession == "ENST00000540868.1"


@pytest.mark.slow
def test_it_does_not_create_entries_for_pseudogenes(human_12):
    entries = {e.optional_id for e in human_12}
    assert "ENSG00000252079.1" not in entries


@pytest.mark.slow
def test_it_normalizes_lineage_to_standard_one(human_12):
    assert entry_for(human_12, "ENST00000540868.1").lineage == (
        "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
        "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; "
        "Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens"
    )


@pytest.mark.slow
def test_calls_lincRNA_lncRNA(human_12):
    assert entry_for(human_12, "ENST00000538041.1").rna_type == "SO:0001877"


@pytest.mark.slow
def test_uses_gene_description_if_possible(human_12):
    assert (
        entry_for(human_12, "ENST00000538041.1").description
        == "Homo sapiens (human) long intergenic non-protein coding RNA 1486"
    )


@pytest.mark.slow
def test_description_strips_source(human_12):
    assert (
        entry_for(human_12, "ENST00000516089.1").description
        == "Homo sapiens (human) small Cajal body-specific RNA 11"
    )


@pytest.mark.slow
def test_generated_description_includes_locus(human_12):
    assert (
        entry_for(human_12, "ENST00000501075.2").description
        == "Homo sapiens (human) novel transcript, antisense to CHD4"
    )


@pytest.mark.slow
def test_can_correct_rfam_name_to_type(human_12):
    assert entry_for(human_12, "ENST00000620330.1").rna_type == "SO:0000590"


@pytest.mark.slow
def test_it_gets_simple_locations(human_12):
    assert entry_for(human_12, "ENST00000546223.1").regions == [
        dat.SequenceRegion(
            chromosome="12",
            strand=-1,
            exons=[
                dat.Exon(start=37773, stop=38102),
                dat.Exon(start=36661, stop=37529),
            ],
            assembly_id="GRCh38",
            coordinate_system=dat.CoordinateSystem.one_based(),
        )
    ]


@pytest.mark.slow
def test_can_get_joined_locations(human_12):
    assert entry_for(human_12, "ENST00000543036.1").regions == [
        dat.SequenceRegion(
            chromosome="12",
            strand=1,
            exons=[
                dat.Exon(start=3319441, stop=3319726),
                dat.Exon(start=3323349, stop=3323452),
                dat.Exon(start=3325090, stop=3325340),
            ],
            assembly_id="GRCh38",
            coordinate_system=dat.CoordinateSystem.one_based(),
        )
    ]


@pytest.mark.slow
def test_it_gets_cross_references(human_12):
    assert entry_for(human_12, "ENST00000504074.1").xref_data == {
        "UCSC": ["ENST00000504074.1"],
        "RNAcentral": ["URS000042090E"],
        "HGNC_trans_name": ["FAM138D-201"],
        "RefSeq_ncRNA": ["NR_026823"],
    }


@pytest.mark.slow
def test_it_uses_correct_antisense_type(human_12):
    assert entry_for(human_12, "ENST00000605233.3").rna_type == "SO:0001877"


@pytest.mark.slow
def test_it_does_not_import_suprressed_rfam_families(human_12):
    assert not entries_for(human_12, "ENST00000611210.1")


@pytest.mark.slow
def test_it_builds_correct_entries(human_12):
    val = attr.asdict(entry_for(human_12, "ENST00000620330.1"))
    del val["sequence"]
    ans = attr.asdict(
        dat.Entry(
            primary_id="ENST00000620330",
            accession="ENST00000620330.1",
            ncbi_tax_id=9606,
            database="ENSEMBL",
            sequence="A",
            regions=[
                dat.SequenceRegion(
                    chromosome="12",
                    strand=1,
                    exons=[dat.Exon(start=3124777, stop=3125063)],
                    assembly_id="GRCh38",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000590",
            url="http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000620330.1",
            seq_version="1",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; "
                "Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens"
            ),
            chromosome="12",
            parent_accession="12.GRCh38",
            common_name="human",
            species="Homo sapiens",
            gene="RF00017",
            locus_tag="RF00017",
            optional_id="ENSG00000278469.1",
            description="Homo sapiens (human) SRP RNA Metazoan signal recognition particle RNA",
            note_data={"transcript_id": ["ENST00000620330.1"]},
            xref_data={
                "UCSC": ["ENST00000620330.1"],
                "RFAM_trans_name": ["RF00017.190-201"],
                "RNAcentral": ["URS0000AA28EF"],
            },
            references=[dat.IdReference(dat.KnownServices.pmid, "27337980")],
            mol_type="genomic DNA",
        )
    )

    del ans["sequence"]
    assert val == ans


@pytest.mark.slow
def test_it_assigns_related(human_x):
    assert entry_for(human_x, "ENST00000434938.7").related_sequences == [
        dat.RelatedSequence(sequence_id="ENST00000430235.7", relationship="isoform")
    ]

    assert entry_for(human_x, "ENST00000430235.7").related_sequences == [
        dat.RelatedSequence(sequence_id="ENST00000434938.7", relationship="isoform")
    ]


@pytest.mark.skip(reason="Not sure is still useful")
def test_it_always_has_valid_rna_types_for_human(human_12):
    for entry in human_12:
        assert entry.rna_type in set(
            [
                "SRP_RNA",
                "Y_RNA",
                "SO:0001904",
                "lncRNA",
                "misc_RNA",
                "other",
                "precursor_RNA",
                "rRNA",
                "ribozyme",
                "snRNA",
                "snoRNA",
                "tRNA",
                "telomerase_RNA",
                "SO:0000584",
            ]
        )


@pytest.mark.slow
def test_it_has_last_ncrna(human_12):
    assert entry_for(human_12, "ENST00000459107.1").xref_data == {
        "RNAcentral": ["URS00006F58F8"],
        "RFAM_trans_name": ["RF00019.633-201"],
        "UCSC": ["ENST00000459107.1"],
    }


@pytest.mark.slow
def test_extracts_all_gencode_entries(human_12):
    assert len([e for e in human_12 if e.database == "GENCODE"]) == 2378


@pytest.mark.slow
def test_can_build_gencode_entries(human_12):
    val = attr.asdict(entry_for(human_12, "GENCODE:ENST00000620330.1"))
    del val["sequence"]
    ans = attr.asdict(
        dat.Entry(
            primary_id="ENST00000620330",
            accession="GENCODE:ENST00000620330.1",
            ncbi_tax_id=9606,
            database="GENCODE",
            sequence="A",
            regions=[
                dat.SequenceRegion(
                    chromosome="12",
                    strand=1,
                    exons=[dat.Exon(start=3124777, stop=3125063)],
                    assembly_id="GRCh38",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000590",
            url="",
            seq_version="1",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; "
                "Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens"
            ),
            chromosome="12",
            parent_accession="12.GRCh38",
            common_name="human",
            species="Homo sapiens",
            gene="RF00017",
            locus_tag="RF00017",
            optional_id=None,
            description="Homo sapiens (human) SRP RNA Metazoan signal recognition particle RNA",
            note_data={"transcript_id": ["ENST00000620330.1"]},
            xref_data={
                "UCSC": ["ENST00000620330.1"],
                "RFAM_trans_name": ["RF00017.190-201"],
                "RNAcentral": ["URS0000AA28EF"],
                "Ensembl": ["ENST00000620330.1"],
            },
            references=[dat.IdReference(dat.KnownServices.pmid, "22955987")],
            mol_type="genomic DNA",
        )
    )

    del ans["sequence"]
    assert val == ans


@pytest.mark.slow
def test_it_does_not_have_excluded_ids(human_12):
    with open("data/ensembl/excluded.txt", "r") as raw:
        excluded = {l.strip() for l in raw}

    # Sanity check
    assert "ENST00000550091.5" in excluded
    assert len(excluded) == 46265

    for entry in human_12:
        name = entry.accession
        if entry.database.lower() == "gencode":
            name = entry.accession.split(":")[1]
        assert entry.accession not in excluded


@pytest.mark.slow
def test_can_use_mouse_models_to_correct_rna_type(mouse_3):
    assert entry_for(mouse_3, "ENSMUST00000082862.1").rna_type == "SO:0000390"


@pytest.mark.skip(reason="Not sure is still useful")
def test_it_always_has_valid_rna_types_for_mouse(mouse_3):
    for entry in mouse_3:
        assert entry.rna_type in set(
            [
                "SRP_RNA",
                "Y_RNA",
                "antisense_RNA",
                "lncRNA",
                "misc_RNA",
                "other",
                "precursor_RNA",
                "rRNA",
                "ribozyme",
                "snRNA",
                "snoRNA",
                "tRNA",
                "SO:0000390",
                "SO:0000584",
            ]
        )


@pytest.mark.slow
def test_correctly_builds_names(human_1):
    assert (
        entry_for(human_1, "ENST00000516935.1").description
        == "Homo sapiens (human) Y RNA"
    )


@pytest.mark.skip(reason="Not sure this is a good test")
def test_it_never_has_bad_vault(mouse_3):
    for entry in mouse_3:
        assert entry.rna_type != "vaultRNA"


@pytest.mark.slow
def test_does_not_append_none_to_description(macaca):
    assert (
        entry_for(macaca, "ENSMMUT00000062476.1").description
        == "Macaca mulatta (rhesus monkey) lncRNA"
    )


@pytest.mark.slow
def test_does_not_create_extra_gencode_entries(macaca):
    assert len([e for e in macaca if e.database == "GENCODE"]) == 0


@pytest.mark.slow
def can_build_reasonable_descriptions_when_locus_is_rfam(cow_8):
    assert (
        entry_for("ENSBTAT00000060095.2").description
        == "Bos taurus (cattle) snoRNA Small nucleolar RNA SNORA8 (ENSBTAG00000043103)"
    )
