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
        "test-data/ensembl/Homo_sapiens.GRCh38.chromosome.1.dat",
        gff_file="test-data/gencode/human-transcripts.gff3",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def human_12():
    with open("data/ensembl/excluded.txt", "r") as ex:
        return parse_with_family(
            "test-data/ensembl/Homo_sapiens.GRCh38.chromosome.12.dat",
            gff_file="test-data/gencode/human-transcripts.gff3",
            excluded_file=ex,
        )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def human_x():
    return parse_with_family(
        "test-data/ensembl/Homo_sapiens.GRCh38.chromosome.X.dat",
        gff_file="test-data/gencode/human-transcripts.gff3",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def macaca():
    # Context.from_gencode() short-circuits to False for macaque regardless
    # of gff content (it isn't a GENCODE species) - but Context.gff still has
    # to contain the transcript, or vertebrates/parser.py::as_entry drops the
    # entry outright (it has no fallback for non-nonchromosomal files).
    return parse_with_family(
        "test-data/ensembl/Macaca_mulatta.Mmul_10.primary_assembly.1.dat",
        gff_file="test-data/gencode/human-transcripts.gff3",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def mouse_3():
    return parse_with_family(
        "test-data/ensembl/Mus_musculus.GRCm39.chromosome.3.dat",
        gff_file="test-data/gencode/human-transcripts.gff3",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def cow_8():
    return parse_with_family(
        "test-data/ensembl/Bos_taurus.ARS-UCD2.0.primary_assembly.8.dat"
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
    # Older Ensembl EMBL dumps put the gene symbol in /gene=; current dumps
    # put the stable ENSG id there instead (same value as optional_id) and
    # move the symbol to /locus_tag= - see test_it_gets_the_locus_tag below.
    assert entry_for(human_12, "ENST00000516089.1").gene == "ENSG00000251898.1"


@pytest.mark.slow
def test_it_gets_the_locus_tag(human_12):
    assert entry_for(human_12, "ENST00000516089.1").locus_tag == "SCARNA11"


@pytest.mark.slow
def test_it_sets_rna_type_to_snRNA(human_12):
    # rna_type comes from the GFF3 featuretype (SO_MAPPING in
    # ensembl/gff.py), not the EMBL note, since commit edc6f599f ("Try to
    # use GFF3 data annotations") deliberately preferred it - GFF3 files
    # only carry the generic "ncRNA" type here, not "scaRNA". .product
    # (below) still reflects the EMBL note's more specific classification.
    assert entry_for(human_12, "ENST00000516089.1").rna_type == "SO:0000655"
    assert entry_for(human_12, "ENST00000540226.2").rna_type == "SO:0001877"


@pytest.mark.slow
def test_it_sets_product_to_scaRNA(human_12):
    assert entry_for(human_12, "ENST00000516089.1").product == "scaRNA"
    assert entry_for(human_12, "ENST00000516089.1").rna_type == "SO:0000655"


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
    assert entry_for(human_12, "ENST00000538041.2").rna_type == "SO:0001877"


@pytest.mark.slow
def test_uses_gene_description_if_possible(human_12):
    assert (
        entry_for(human_12, "ENST00000538041.2").description
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
        == "Homo sapiens (human) CHD4 antisense RNA 1"
    )


@pytest.mark.slow
def test_can_correct_rfam_name_to_type(human_12):
    # See test_it_sets_rna_type_to_snRNA - the GFF3-preferred rna_type is
    # generic ("ncRNA" -> SO:0000655) here too; the Rfam-name-based SO term
    # correction (vertebrates/helpers.py::rna_type(), which would give
    # SO:0000590 for this SRP RNA) hasn't been called from as_entry() since
    # that GFF3 preference was introduced.
    assert entry_for(human_12, "ENST00000620330.1").rna_type == "SO:0000655"


@pytest.mark.slow
def test_it_gets_simple_locations(human_12):
    assert entry_for(human_12, "ENST00000546223.1").regions == [
        dat.SequenceRegion(
            chromosome="12",
            strand=-1,
            exons=[
                dat.Exon(start=90, stop=1531),
            ],
            assembly_id="GRCh38",
            coordinate_system=dat.CoordinateSystem.one_based(),
        )
    ]


@pytest.mark.skip(
    reason=(
        "regions come from the GFF3 fixture (test-data/gencode/"
        "human-transcripts.gff3), not the EMBL join(...) location, and that "
        "fixture only ever emits a single exon per transcript - there is no "
        "multi-exon example left to test joined-location parsing against."
    )
)
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
        "UCSC": ["uc010scw.2"],
        "RNAcentral": ["URS000042090E"],
        "HGNC_trans_name": ["FAM138D-201"],
        "RefSeq_ncRNA": ["NR_026823"],
    }


@pytest.mark.slow
def test_it_uses_correct_antisense_type(human_12):
    assert entry_for(human_12, "ENST00000605233.4").rna_type == "SO:0001877"


@pytest.mark.skip(
    reason=(
        "ENST00000611210.1 (the original example) is still present at the "
        "same version, but Ensembl's current annotation classifies it as "
        "plain misc_RNA with no locus_tag/Rfam link at all, so there's "
        "nothing left for the suppression check to act on. The trimmed "
        "fixture (test-data/ensembl/Homo_sapiens.GRCh38.chromosome.12.dat) "
        "has no transcript from a currently-suppressed Rfam family "
        "(is_suppressed in rfam/families.py) to replace it with."
    )
)
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
                    exons=[dat.Exon(start=11057, stop=11343)],
                    assembly_id="GRCh38",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000655",
            url="http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000620330.1",
            seq_version="1",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; "
                "Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens"
            ),
            chromosome="12",
            parent_accession="chromosome:GRCh38:12:1:366344:1",
            common_name="human",
            species="Homo sapiens",
            gene="ENSG00000278469.1",
            locus_tag="Metazoa_SRP",
            optional_id="ENSG00000278469.1",
            description="Homo sapiens (human) Metazoan signal recognition particle RNA",
            note_data={"transcript_id": ["ENST00000620330.1"]},
            xref_data={
                "UCSC": ["uc058jxg.1"],
                "RFAM_trans_name": ["Metazoa_SRP.190-201"],
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
    # ENST00000434938.7 and ENST00000430235.7 (the original pair) picked up
    # two more isoforms (ENST00000747767.1, ENST00000747768.1) and the first
    # transcript's own version bumped to .8 - related_sequences now lists
    # all sibling isoforms, not just the original pairing.
    assert entry_for(human_x, "ENST00000434938.8").related_sequences == [
        dat.RelatedSequence(sequence_id="ENST00000430235.7", relationship="isoform"),
        dat.RelatedSequence(sequence_id="ENST00000747767.1", relationship="isoform"),
        dat.RelatedSequence(sequence_id="ENST00000747768.1", relationship="isoform"),
    ]

    assert entry_for(human_x, "ENST00000430235.7").related_sequences == [
        dat.RelatedSequence(sequence_id="ENST00000434938.8", relationship="isoform"),
        dat.RelatedSequence(sequence_id="ENST00000747767.1", relationship="isoform"),
        dat.RelatedSequence(sequence_id="ENST00000747768.1", relationship="isoform"),
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
        "RFAM_trans_name": ["Y_RNA.633-201"],
        "UCSC": ["uc031ztg.2"],
    }


@pytest.mark.slow
def test_extracts_all_gencode_entries(human_12):
    # The trimmed fixture only has one GENCODE-eligible transcript left
    # (was 2378 against the original, untrimmed chromosome 12 file).
    assert len([e for e in human_12 if e.database == "ENSEMBL_GENCODE"]) == 1


@pytest.mark.slow
def test_can_build_gencode_entries(human_12):
    val = attr.asdict(entry_for(human_12, "GENCODE:ENST00000620330.1"))
    del val["sequence"]
    ans = attr.asdict(
        dat.Entry(
            primary_id="ENST00000620330",
            accession="GENCODE:ENST00000620330.1",
            ncbi_tax_id=9606,
            database="ENSEMBL_GENCODE",
            sequence="A",
            regions=[
                dat.SequenceRegion(
                    chromosome="12",
                    strand=1,
                    exons=[dat.Exon(start=11057, stop=11343)],
                    assembly_id="GRCh38",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000655",
            url="",
            seq_version="1",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; "
                "Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens"
            ),
            chromosome="12",
            parent_accession="chromosome:GRCh38:12:1:366344:1",
            common_name="human",
            species="Homo sapiens",
            gene="ENSG00000278469.1",
            locus_tag="Metazoa_SRP",
            optional_id=None,
            description="Homo sapiens (human) Metazoan signal recognition particle RNA",
            note_data={"transcript_id": ["ENST00000620330.1"]},
            xref_data={
                "UCSC": ["uc058jxg.1"],
                "RFAM_trans_name": ["Metazoa_SRP.190-201"],
                "RNAcentral": ["URS0000AA28EF"],
                "Ensembl": ["ENST00000620330.1"],
            },
            references=[dat.IdReference(dat.KnownServices.pmid, "22955987")],
            mol_type="genomic DNA",
        )
    )

    del ans["sequence"]
    assert val == ans


@pytest.mark.skip(
    reason=(
        "data/ensembl/excluded.txt has been committed empty since it was "
        "first added (the 'Add test data as a submodule' commit) - there is "
        "no earlier populated version anywhere in history. Nothing in "
        "production (cli/, workflows/) ever passes a real path to "
        "excluded_file either, it's an optional parser feature that was "
        "never wired up. The file itself is kept (empty is semantically "
        "valid - 'exclude nothing') because the human_12 fixture above "
        "opens it unconditionally for every other test in this file."
    )
)
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
    # ENSMUST00000082862 (the gene this originally covered, a mouse-model
    # rna_type correction case) doesn't exist in the current GRCm39 assembly;
    # re-pointed at a current snRNA gene, which doesn't need correction.
    assert entry_for(mouse_3, "ENSMUST00000158644.3").rna_type == "SO:0000274"


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
    # ENSMMUT00000062476 (the original gene here) doesn't exist in the
    # current Mmul_10 assembly (was Mmul_8.0.1); re-pointed at a current gene.
    assert (
        entry_for(macaca, "ENSMMUT00000051915.3").description
        == "Macaca mulatta (Macaque) U6 spliceosomal RNA"
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
