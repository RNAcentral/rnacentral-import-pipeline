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

from io import open
from pathlib import Path

import attr
import pytest

import rnacentral_pipeline.databases.helpers.publications as pubs
from rnacentral_pipeline.databases.data import Entry, Reference
from rnacentral_pipeline.databases.ena import context, parser, ribovore


def simple_parse(path):
    builder = context.ContextBuilder()
    builder.with_dr(path)
    ctx = builder.context()
    return parser.parse(ctx, path)


@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/ena/ncr/wgs/aa/wgs_aacd01_fun.ncr", 188),
        ("data/ena/ncr/wgs/aa/wgs_abxv02_pro.ncr", 103),
        ("data/test_AS_lines.ncr", 1),
        ("data/test_example_entries.ncr", 10),
        ("data/test_feature_unclosed_parenthesis.ncr", 2),
        ("data/test_invalid_fields.ncr", 1),
        ("data/test_invalid_sequences.ncr", 4),
        ("data/test_long_sequence.ncr", 1),
        ("data/test_rdp_product.ncr", 1),
        ("data/test_rfam_duplicates.ncr", 8),
        ("data/test_rfam_product.ncr", 2),  # Skipped because outdated taxon id
        ("data/test_silva_product.ncr", 3),
        ("data/test_species_patch.ncr", 1),
        ("data/ena/ncr/ex/wgs_acnt01_pro.ncr", 105),
        ("data/ena/protein-to-skip.embl", 0),  # Proteins that should be skipped
        ("data/ena/mislabeled-trna.embl", 1),
        ("data/ena/too-long-description.ncr", 1),
        ("data/ena/scarna.ncr", 3),
        ("data/ena/bad-feature-type.ncr", 2),
    ],
)
def test_can_parse_variety_of_files(filename, count):
    data = list(simple_parse(Path(filename)))
    assert len(data) == count
    for entry in data:
        assert "TPA:" not in entry.description
        assert entry.description
        for reference in entry.references:
            if not isinstance(reference, Reference):
                continue
            assert reference.title != ";"


def test_creates_simple_entry():
    path = Path("data/ena/ncr/wgs/aa/wgs_aacd01_fun.ncr")
    simple_data = list(simple_parse(path))

    assert attr.asdict(simple_data[0]) == attr.asdict(
        Entry(
            primary_id="",
            accession="AACD01000002.1:101667..101773:tRNA",
            ncbi_tax_id=227321,
            database="ENA",
            sequence="GCCCGGATGGTCTAGTGGTATGATTCTCCCTTCGGGGGCAGCGCCCGGTACATAATAACATGTATCAGAAATGGGAGAGGTCCCGCGTTCGAATCGCGGTTCGGGCC",
            regions=[],
            rna_type="SO:0000253",  # Ideally SO:0000268
            url="https://www.ebi.ac.uk/ena/data/view/Non-coding:AACD01000002.1:101667..101773:tRNA",
            seq_version="1",
            xref_data={"ena_refs": {"BIOSAMPLE": ("SAMN02953587", None)}},
            chromosome="VIII",
            species="Aspergillus nidulans FGSC A4",
            lineage=(
                "Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; "
                "Eurotiomycetes; Eurotiomycetidae; Eurotiales; Aspergillaceae; "
                "Aspergillus; Aspergillus nidulans FGSC A4"
            ),
            product="tRNA-Pro",
            parent_accession="AACD01000002",
            description="Aspergillus nidulans FGSC A4 tRNA-Pro",
            mol_type="genomic DNA",
            references=[
                pubs.reference(16372000),
                Reference(
                    location=(
                        "Submitted (04-APR-2003) to the INSDC. Whitehead "
                        "Institute/MIT Center for Genome Research, 320 Charles "
                        "Street, Cambridge, MA 02142, USA"
                    ),
                    title=None,
                    authors=(
                        "Birren B., Nusbaum C., Abouelleil A., Allen N., "
                        "Anderson S., Arachchi H.M., Barna N., Bastien V., "
                        "Bloom T., Boguslavkiy L., Boukhgalter B., Butler J., "
                        "Calvo S.E., Camarata J., Chang J., Choepel Y., "
                        "Collymore A., Cook A., Cooke P., Corum B., DeArellano K.,"
                        " Diaz J.S., Dodge S., Dooley K., Dorris L., Elkins T., "
                        "Engels R., Erickson J., Faro S., Ferreira P., "
                        "FitzGerald M., Gage D., Galagan J., Gardyna S., Gnerre S."
                        ", Graham L., Grand-Pierre N., Hafez N., Hagopian D., "
                        "Hagos B., Hall J., Horton L., Hulme W., Iliev I., "
                        "Jaffe D., Johnson R., Jones C., Kamal M., Kamat A., "
                        "Karatas A., Kells C., Landers T., Levine R., "
                        "Lindblad-Toh K., Liu G., Lui A., Ma L.-J., Mabbitt R., "
                        "MacLean C., Macdonald P., Major J., Manning J., "
                        "Matthews C., Mauceli E., McCarthy M., Meldrim J., "
                        "Meneus L., Mihova T., Mlenga V., Murphy T., Naylor J., "
                        "Nguyen C., Nicol R., Nielsen C.B., Norbu C., O'Connor T.,"
                        " O'Donnell P., O'Neil D., Oliver J., Peterson K., "
                        "Phunkhang P., Pierre N., Purcell S., Rachupka A., "
                        "Ramasamy U., Raymond C., Retta R., Rise C., Rogov P., "
                        "Roman J., Schauer S., Schupback R., Seaman S., Severy P.,"
                        " Smirnov S., Smith C., Spencer B., Stange-Thomann N., "
                        "Stojanovic N., Stubbs M., Talamas J., Tesfaye S., "
                        "Theodore J., Topham K., Travers M., Vassiliev H., "
                        "Venkataraman V.S., Viel R., Vo A., Wang S., Wilson B., "
                        "Wu X., Wyman D., Young G., Zainoun J., Zembek L., "
                        "Zimmer A., Zody M., Lander E."
                    ),
                    pmid=None,
                    doi=None,
                ),
                Reference(
                    authors=(
                        "Birren B., Nusbaum C., Abebe A., Abouelleil A., Adekoya "
                        "E., Ait-zahra M., Allen N., Allen T., An P., Anderson M., "
                        "Anderson S., Arachchi H., Armbruster J., Bachantsang P., "
                        "Baldwin J., Barry A., Bayul T., Blitshsteyn B., Bloom T., "
                        "Blye J., Boguslavskiy L., Borowsky M., Boukhgalter B., "
                        "Brunache A., Butler J., Calixte N., Calvo S., Camarata "
                        "J., Campo K., Chang J., Cheshatsang Y., Citroen M., "
                        "Collymore A., Considine T., Cook A., Cooke P., Corum B., "
                        "Cuomo C., David R., Dawoe T., Degray S., Dodge S., Dooley "
                        "K., Dorje P., Dorjee K., Dorris L., Duffey N., Dupes A., "
                        "Elkins T., Engels R., Erickson J., Farina A., Faro S., "
                        "Ferreira P., Fischer H., Fitzgerald M., Foley K., Gage "
                        "D., Galagan J., Gearin G., Gnerre S., Gnirke A., Goyette "
                        "A., Graham J., Grandbois E., Gyaltsen K., Hafez N., "
                        "Hagopian D., Hagos B., Hall J., Hatcher B., Heller A., "
                        "Higgins H., Honan T., Horn A., Houde N., Hughes L., Hulme "
                        "W., Husby E., Iliev I., Jaffe D., Jones C., Kamal M., "
                        "Kamat A., Kamvysselis M., Karlsson E., Kells C., Kieu A., "
                        "Kisner P., Kodira C., Kulbokas E., Labutti K., Lama D., "
                        "Landers T., Leger J., Levine S., Lewis D., Lewis T., "
                        "Lindblad-toh K., Liu X., Lokyitsang T., Lokyitsang Y., "
                        "Lucien O., Lui A., Ma L.J., Mabbitt R., Macdonald J., "
                        "Maclean C., Major J., Manning J., Marabella R., Maru K., "
                        "Matthews C., Mauceli E., Mccarthy M., Mcdonough S., "
                        "Mcghee T., Meldrim J., Meneus L., Mesirov J., Mihalev A., "
                        "Mihova T., Mikkelsen T., Mlenga V., Moru K., Mozes J., "
                        "Mulrain L., Munson G., Naylor J., Newes C., Nguyen C., "
                        "Nguyen N., Nguyen T., Nicol R., Nielsen C., Nizzari M., "
                        "Norbu C., Norbu N., O'donnell P., Okoawo O., O'leary S., "
                        "Omotosho B., O'neill K., Osman S., Parker S., Perrin D., "
                        "Phunkhang P., Piqani B., Purcell S., Rachupka T., "
                        "Ramasamy U., Rameau R., Ray V., Raymond C., Retta R., "
                        "Richardson S., Rise C., Rodriguez J., Rogers J., Rogov "
                        "P., Rutman M., Schupbach R., Seaman C., Settipalli S., "
                        "Sharpe T., Sheridan J., Sherpa N., Shi J., Smirnov S., "
                        "Smith C., Sougnez C., Spencer B., Stalker J., "
                        "Stange-thomann N., Stavropoulos S., Stetson K., Stone C., "
                        "Stone S., Stubbs M., Talamas J., Tchuinga P., Tenzing P., "
                        "Tesfaye S., Theodore J., Thoulutsang Y., Topham K., Towey "
                        "S., Tsamla T., Tsomo N., Vallee D., Vassiliev H., "
                        "Venkataraman V., Vinson J., Vo A., Wade C., Wang S., "
                        "Wangchuk T., Wangdi T., Whittaker C., Wilkinson J., Wu "
                        "Y., Wyman D., Yadav S., Yang S., Yang X., Yeager S., Yee "
                        "E., Young G., Zainoun J., Zembeck L., Zimmer A., Zody M., "
                        "Lander E."
                    ),
                    location=(
                        "Submitted (07-JAN-2004) to the INSDC. Whitehead "
                        "Institute/MIT Center for Genome Research, "
                        "320 Charles Street, Cambridge, MA 02142, USA"
                    ),
                    title=None,
                    pmid=None,
                    doi=None,
                ),
            ],
        )
    )


def test_can_find_correct_ncRNA_type():
    path = Path("data/ena/ncr/wgs/aa/wgs_abxv02_pro.ncr")
    ncrna_data = list(simple_parse(path))

    assert attr.asdict(ncrna_data[0]) == attr.asdict(
        Entry(
            primary_id="",
            accession="ABXV02000002.1:33498..33573:ncRNA",
            ncbi_tax_id=500637,
            database="ENA",
            sequence="ACTGCTTTTCTTTGATGTCCCCATATTGAGGAGCCCGATAGCCATTTGATTACTTCATGCTATCGGGTTTTTTATT",
            regions=[],
            rna_type="SO:0000655",
            url="https://www.ebi.ac.uk/ena/data/view/Non-coding:ABXV02000002.1:33498..33573:ncRNA",
            seq_version="1",
            note_data={"text": ["Rfam score 66"]},
            xref_data={"ena_refs": {"BIOSAMPLE": ("SAMN00008855", None)}},
            species="Providencia rustigianii DSM 4541",
            lineage=(
                "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; "
                "Morganellaceae; Providencia; Providencia rustigianii DSM 4541"
            ),
            product="RybB RNA",
            parent_accession="ABXV02000002",
            description="Providencia rustigianii DSM 4541 RybB RNA",
            mol_type="genomic DNA",
            inference="nucleotide motif:Rfam:RF00110",
            locus_tag="PROVRUST_04548",
            references=[
                Reference(
                    authors=(
                        "Fulton L., Clifton S., Fulton B., Xu J., Minx P., Pepin "
                        "K.H., Johnson M., Bhonagiri V., Nash W.E., Mardis E.R., "
                        "Wilson R.K."
                    ),
                    location=(
                        "Submitted (16-OCT-2008) to the INSDC. Genome Sequencing "
                        "Center, Washington University School of Medicine, 4444 "
                        "Forest Park, St. Louis, MO 63108, USA"
                    ),
                    title=None,
                    pmid=None,
                    doi=None,
                ),
                Reference(
                    authors=(
                        "Weinstock G., Sodergren E., Clifton S., Fulton L., Fulton "
                        "B., Courtney L., Fronick C., Harrison M., Strong C., "
                        "Farmer C., Delahaunty K., Markovic C., Hall O., Minx P., "
                        "Tomlinson C., Mitreva M., Nelson J., Hou S., Wollam A., "
                        "Pepin K.H., Johnson M., Bhonagiri V., Nash W.E., Warren "
                        "W., Chinwalla A., Mardis E.R., Wilson R.K."
                    ),
                    location=(
                        "Submitted (21-SEP-2009) to the INSDC. Genome Sequencing "
                        "Center, Washington University School of Medicine, 4444 "
                        "Forest Park, St. Louis, MO 63108, USA"
                    ),
                    title=None,
                    pmid=None,
                    doi=None,
                ),
                Reference(
                    authors=(
                        "Weinstock G., Sodergren E., Clifton S., Fulton L., Fulton "
                        "B., Courtney L., Fronick C., Harrison M., Strong C., "
                        "Farmer C., Delahaunty K., Markovic C., Hall O., Minx P., "
                        "Tomlinson C., Mitreva M., Nelson J., Hou S., Wollam A., "
                        "Pepin K.H., Johnson M., Bhonagiri V., Nash W.E., Warren "
                        "W., Chinwalla A., Mardis E.R., Wilson R.K."
                    ),
                    location=(
                        "Submitted (14-DEC-2009) to the INSDC. Genome Sequencing "
                        "Center, Washington University School of Medicine, 4444 "
                        "Forest Park, St. Louis, MO 63108, USA"
                    ),
                    title=None,
                    pmid=None,
                    doi=None,
                ),
            ],
        )
    )


def test_can_parse_all_example_entries():
    path = Path("data/test_example_entries.ncr")
    examples = list(simple_parse(path))

    assert attr.asdict(examples[0]) == attr.asdict(
        Entry(
            primary_id="",
            accession="AB330787.1:1..34:misc_RNA",
            ncbi_tax_id=9606,
            database="ENA",
            sequence="ATTGGGGAGTGAGAAGGAGAGAACGCGGTCTGAA",
            regions=[],
            rna_type="SO:0000673",
            url="https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330787.1:1..34:misc_RNA",
            seq_version="1",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
                "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
                "Homo sapiens"
            ),
            species="Homo sapiens",
            common_name="human",
            description="Homo sapiens (human) miscellaneous RNA",
            parent_accession="AB330787",
            note_data={
                "text": [
                    "RNA sequence binds to tRNaseZL",
                    "similar to U3 snRNA fragment",
                ]
            },
            mol_type="other RNA",
        )
    )

    assert attr.asdict(examples[1]) == attr.asdict(
        Entry(
            primary_id="",
            accession="AB330786.1:1..27:misc_RNA",
            database="ENA",
            ncbi_tax_id=9606,
            sequence="ATTGCAGTACCTCCAGGAACGGTGCAC",
            regions=[],
            rna_type="SO:0000673",
            url="https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330786.1:1..27:misc_RNA",
            seq_version="1",
            xref_data={"ena_refs": {"RFAM": ("RF00005", "tRNA")}},
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
                "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
                "Homo sapiens"
            ),
            species="Homo sapiens",
            common_name="human",
            description="Homo sapiens (human) miscellaneous RNA",
            parent_accession="AB330786",
            note_data={
                "text": [
                    "RNA sequence binds to tRNaseZL",
                    "similar to U2 snRNA fragment",
                ]
            },
            mol_type="other RNA",
        )
    )

    assert attr.asdict(examples[2]) == attr.asdict(
        Entry(
            primary_id="",
            accession="AB330785.1:1..34:misc_RNA",
            database="ENA",
            ncbi_tax_id=9606,
            sequence="CGCGACCTCAGATCAGACGTGGCGACCCGCTGAA",
            regions=[],
            rna_type="SO:0000673",
            url="https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330785.1:1..34:misc_RNA",
            seq_version="1",
            xref_data={"ena_refs": {"SRPDB": ("Xeno.laev._DC015625", None)}},
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
                "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
                "Homo sapiens"
            ),
            species="Homo sapiens",
            common_name="human",
            description="Homo sapiens (human) miscellaneous RNA",
            note_data={
                "text": [
                    "RNA sequence binds to tRNaseZL",
                    "similar to 28S rRNA fragment",
                ]
            },
            parent_accession="AB330785",
            mol_type="other RNA",
        )
    )

    assert attr.asdict(examples[3]) == attr.asdict(
        Entry(
            primary_id="",
            accession="HAAO01001079.1:1..21:ncRNA",
            database="ENA",
            ncbi_tax_id=9606,
            sequence="ACCACTGCACTCCAGCCTGAG",
            regions=[],
            rna_type="miRNA",
            url="https://www.ebi.ac.uk/ena/data/view/Non-coding:HAAO01001079.1:1..21:ncRNA",
            seq_version="1",
            xref_data={"ena_refs": {"MIRBASE": ("MIMAT0022742", "hsa-miR-1273g-3p")}},
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
                "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
                "Homo sapiens"
            ),
            species="Homo sapiens",
            common_name="human",
            description="Homo sapiens (human) microRNA hsa-miR-1273g-3p",
            product="microRNA hsa-miR-1273g-3p",
            mol_type="transcribed RNA",
            gene="hsa-miR-1273g-3p",
            parent_accession="HAAO01001079",
            references=[
                Reference(
                    authors="",
                    location="Submitted (20-AUG-2013) to the INSDC.",
                    title=None,
                    pmid=None,
                    doi=None,
                ),
                pubs.reference(21037258),
            ],
        )
    )

    assert attr.asdict(examples[4]) == attr.asdict(
        Entry(
            primary_id="",
            accession="HG519048.1:1..359:tmRNA",
            database="ENA",
            ncbi_tax_id=32644,
            sequence="GGGAGCGACTTGGCTTCGACAGGAGTAAGTCTGCTTAGATGGCATGTCGCTTTGGGCAAAGCGTAAAAAGCCCAAATAAAATTAAACGCAAACAACGTTAAATTCGCTCCTGCTTACGCTAAAGCTGCGTAAGTTCAGTTGAGCCTGAAATTTAAGTCATACTATCTAGCTTAATTTTCGGTCATTTTTGATAGTGTAGCCTTGCGTTTGACAAGCGTTGAGGTGAAATAAAGTCTTAGCCTTGCTTTTGAGTTTTGGAAGATGAGCGAAGTAGGGTGAAGTAGTCATCTTTGCTAAGCATGTAGAGGTCTTTGTGGGATTATTTTTGGACAGGGGTTCGATTCCCCTCGCTTCCACCA",
            regions=[],
            rna_type="tmRNA",
            url="https://www.ebi.ac.uk/ena/data/view/Non-coding:HG519048.1:1..359:tmRNA",
            seq_version="1",
            note_data={
                "text": ["Tag:(A)ANNVKFAPAYAKAA*"],
            },
            xref_data={"ena_refs": {"TMRNA-WEBSITE": ("Campy_jejun_700819", None)}},
            lineage="unclassified sequences; unidentified",
            species="unidentified",
            description="unidentified transfer-messenger mRNA Campy_jejun_700819",
            product="transfer-messenger mRNA Campy_jejun_700819",
            gene="tmRNA Campy_jejun_700819",
            mol_type="genomic DNA",
            parent_accession="HG519048",
            references=[
                Reference(
                    authors="",
                    location="Submitted (13-SEP-2013) to the INSDC.",
                    title=None,
                    pmid=None,
                    doi=None,
                ),
                pubs.reference(14681369),
            ],
        )
    )

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='HG975377.1:1..332:ncRNA',
    #     database='ENA',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='HG985290.1:1..72:tRNA',
    #     database='ENA',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='LM608264.1:7..26:ncRNA',
    #     database='ENA',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='BX284601.5:1693190..1693457:ncRNA',
    #     database='ENA',
    #     rna_type='SO:0000275',
    #     product="RNA transcript Y71G12B.41",
    #     standard_name="Y71G12B.41",
    #     locus_tag="CELE_Y71G12B.41",
    #     gene="Y71G12B.41",
    #     mol_type='genomic DNA',
    #     chromosome='I',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='CU329672.1:1601457..1602924:misc_RNA',
    #     database='ENA',
    #     ncbi_tax_id=4896,
    #     rna_type='SO:0000673'
    #     mol_type='genomic DNA',
    #     locus_tag="SPNCRNA.1210",
    #     chromosome='III',
    #     product='antisense RNA (predicted)'
    # )),


def test_can_handle_file_with_invalid_fields():
    raw = Path("data/test_invalid_fields.ncr")
    data = next(simple_parse(raw))

    assert data.accession == "KM079256.1:1..1300:rRNA"
    assert data.mol_type == "genomic DNA"
    assert data.species == 'Candidatus Stammerula sp. of Trupanea "pohakuloa"'
    assert data.ncbi_tax_id == 1630665
    assert data.product == "16S ribosomal RNA"


@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/ena/pseudogene.embl", 0),
        ("data/ena/pseudogene-flag.embl", 0),
    ],
)
def correctly_ignores_pseudogenes(filename, count):
    raw = Path(filename)
    data = list(simple_parse(raw))
    assert len(data) == count
    # assert data.accession == 'NIDN01000248.1:29758..29889:tRNA'
    # assert data.pseudogene == 'unprocessed'


def test_can_parse_function():
    raw = Path("data/ena/function.embl")
    data = list(simple_parse(raw))

    assert data[0].accession == "EU410654.1:1..92:ncRNA"
    assert data[0].function == "guide for 26S rRNA methylation at U1043"
    assert data[0].references[0] == pubs.reference(18493037)
    assert data[1].accession == "AB046489.1:221..306:tRNA"
    assert data[1].function == "tRNA-Pro"
    assert data[1].organelle == "mitochondrion"

    assert data[2].accession == "CP003783.1:1548698..1548818:ncRNA"
    assert data[2].function == "1.8: Sporulation"
    assert data[2].gene == "csfG"
    assert data[2].product == "sporulation-specific regulatory RNA"
    assert data[2].locus_tag == "B657_miscRNA23"
    assert data[2].rna_type == "SO:0000655"


def test_can_parse_gene_synonyms():
    raw = Path("data/ena/gene_synonym.embl")
    data = next(simple_parse(raw))

    assert data.accession == "CP000948.1:2011661..2011909:misc_RNA"
    assert data.gene_synonyms == ["IS091", "sraC", "tpke79"]
    assert data.locus_tag == "ECDH10B_1978"
    assert data.mol_type == "genomic DNA"
    assert data.gene == "ryeA"
    assert data.product == "small RNA"


@pytest.mark.parametrize(
    "pmid",
    [
        12244299,
        6196367,
        6181418,
        6802847,
        3403542,
        6084597,
        10924331,
        7528809,
        7529207,
        1704372,
        20610725,
        18617187,
        17881443,
        17164479,
        8389475,
        10834842,
        10684931,
        15611297,
        20668672,
        911771,
        6209580,
    ],
)
def test_can_extract_references_from_experiment(pmid):
    raw = Path("data/ena/experiment-references.embl")
    data = next(simple_parse(raw))

    assert data.accession == "HG975378.1:1..299:ncRNA"
    assert data.rna_type == "SO:0000590"
    assert data.gene == "RN7SL1"
    assert data.mol_type == "transcribed RNA"
    assert data.product == "Small nucleolar RNA 7SL"
    assert data.note_data == {
        "ontology": ["ECO:0000305", "GO:0006617", "GO:0048501", "SO:0000590"],
        "text": ["biotype:SRP_RNA"],
    }
    assert pubs.reference(pmid) in data.references


@pytest.mark.parametrize(
    "pmid", [1379177, 8288542, 9098041, 9826762, 12165569, 12547201, 1317842, 15831787]
)
def test_can_extract_references_from_note_or_experiment(pmid):
    raw = Path("data/ena/note-references.embl")
    data = next(simple_parse(raw))

    assert data.accession == "AL009126.3:3856226..3856447:misc_RNA"
    assert data.note_data == {
        "text": [
            "Evidence 1a: Function experimentally demonstrated in the studied strain; PubMedId: 1379177, 8288542, 9098041, 9826762, 12165569, 12547201, 1317842, 15831787; Product type n: RNA"
        ]
    }
    assert data.product == "T-box"
    assert data.locus_tag == "BSU_misc_RNA_57"
    assert data.function == "16.3: Control"
    assert data.gene == "tboTB"
    assert pubs.reference(pmid) in data.references


def test_can_handle_unclosed_parens():
    raw = Path("data/test_feature_unclosed_parenthesis.ncr")
    data = next(simple_parse(raw))

    assert data.accession == "HE860504.1:1..14644:tRNA"
    assert data.mol_type == "genomic DNA"
    assert data.ncbi_tax_id == 1200666
    assert data.gene == "tRNA-Ser"
    assert data.product == "transfer RNA Serine"
    assert data.species == "Metacrangonyx sp. 3 ssp. 1 MDMBR-2012"
    assert data.organelle == "mitochondrion"
    assert (
        data.description == "Metacrangonyx sp. 3 ssp. 1 MDMBR-2012 transfer RNA Serine"
    )


def test_can_extract_xrefs():
    raw = Path("data/test_example_entries.ncr")
    data = list(simple_parse(raw))

    assert data[5].accession == "HG975377.1:1..332:ncRNA"
    assert data[5].xref_data == {"lncrnadb": ["190", "7SK"]}
    assert data[7].accession == "LM608264.1:7..26:ncRNA"
    assert data[7].xref_data == {"mirbase": ["MI0000182"]}


def test_can_handle_mislabled_trna():
    raw = Path("data/ena/mislabeled-trna.embl")
    data = list(simple_parse(raw))

    assert data[0].rna_type == "SO:0000253"


@pytest.mark.parametrize(
    "index,rna_type",
    [
        (0, "SO:0000252"),
        (1, "SO:0000252"),
    ],
)
def test_can_handle_weird_feature(index, rna_type):
    raw = Path("data/ena/bad-feature-type.ncr")
    data = list(simple_parse(raw))

    assert data[index].rna_type == rna_type


@pytest.mark.parametrize(
    "index,rna_type",
    [
        (0, "SO:0000602"),
        (1, "SO:0000275"),
        (2, "SO:0000602"),
    ],
)
def test_can_get_scarna_feature(index, rna_type):
    raw = Path("data/ena/scarna.ncr")
    data = list(simple_parse(raw))

    assert data[index].rna_type == rna_type


def test_deals_with_crazy_long_name():
    raw = Path("data/ena/too-long-description.ncr")
    data = list(simple_parse(raw))

    assert data[0].rna_type == "SO:0000275"
    assert data[0].description == "Leishmania donovani snoRNA.05.001"
    assert data[0].locus_tag == "LDHU3_LDHU3_05.T1400"
    # assert data[0].product == 'snoRNA'
    # assert data[0].gene == 'snoRNA.05.001'
    # assert data[0].gene_synonyms == [
    #     'snoRNA.05.002',
    #     'snoRNA.05.003',
    #     'snoRNA.05.004',
    #     'snoRNA.05.005',
    #     'snoRNA.05.006',
    #     'snoRNA.05.007',
    #     'snoRNA.05.008',
    #     'snoRNA.05.009',
    #     'snoRNA.05.010',
    #     'snoRNA.05.011',
    #     'snoRNA.05.012',
    #     'snoRNA.05.013',
    #     'snoRNA.05.014',
    #     'snoRNA.05.014',
    #     'snoRNA.05.015',
    #     'snoRNA.05.016',
    #     'snoRNA.05.017',
    #     'snoRNA.05.018',
    #     'snoRNA.05.019',
    #     'snoRNA.05.020',
    #     'snoRNA.05.021',
    #     'snoRNA.05.022',
    #     'snoRNA.05.023',
    #     'snoRNA.05.024',
    #     'snoRNA.05.025',
    #     'snoRNA.05.026',
    #     'snoRNA.05.027',
    # ]


@pytest.mark.parametrize(
    "path,id,included",
    [
        ("data/ena/to-exclude/lazr01007255", "LAZR01007255.1:2..3402:rRNA", True),
        ("data/ena/to-exclude/snry01001410", "SNRY01001410.1:2783..5704:rRNA", True),
        pytest.param(
            "data/ena/to-exclude/haqp01000579",
            "HAQP01000579.1:18..191:rRNA",
            False,
            marks=pytest.mark.xfail(reason="TBI"),
        ),
        ("data/ena/to-exclude/kbtv01008406", "KBTV01008406.1:1..998:rRNA", False),
        ("data/ena/to-exclude/kdhv01000321", "KDHV01000321.1:1..277:rRNA", False),
        ("data/ena/to-exclude/kdvk01017574", "KDVK01017574.1:1..441:rRNA", False),
    ],
)
def test_knows_how_to_exclude_or_not(path, id, included):
    path = Path(path)
    raw = path / "sequences.dat"
    ribotyper = path / "ribotyper-results"
    model_lengths = Path("data/ena/to-exclude/model-lengths.csv")
    builder = context.ContextBuilder()
    builder.with_dr(raw)
    builder.with_ribovore(ribotyper, model_lengths)
    ctx = builder.context()
    data = list(parser.parse(ctx, raw))
    was_emitted = any(e.accession == id for e in data)
    assert was_emitted == included
