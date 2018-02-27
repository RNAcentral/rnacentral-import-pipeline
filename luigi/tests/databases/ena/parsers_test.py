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

from databases.data import Entry, Exon, Reference
from databases.ena.parsers import parse


@pytest.mark.parametrize('filename,count', [
    ('data/ena/ncr/wgs/aa/wgs_aacd01_fun.ncr', 188),
    ('data/ena/ncr/wgs/aa/wgs_abxv02_pro.ncr', 103),
    ('data/test_AS_lines.ncr', 1),
    ('data/test_example_entries.ncr', 10),
    ('data/test_feature_unclosed_parenthesis.ncr', 2),
    ('data/test_invalid_fields.ncr', 1),
    ('data/test_invalid_sequences.ncr', 4),
    ('data/test_long_sequence.ncr', 1),
    ('data/test_rdp_product.ncr', 1),
    ('data/test_refseq_product.ncr', 3),
    ('data/test_rfam_duplicates.ncr', 8),
    ('data/test_rfam_product.ncr', 2),  # Skipped because outdated taxon id
    ('data/test_silva_product.ncr', 3),
    ('data/test_species_patch.ncr', 1),
    ('data/ena/ncr/ex/wgs_acnt01_pro.ncr', 106),
])
def test_can_parse_variety_of_files(filename, count):
    with open(filename, 'rb') as raw:
        print(filename)
        data = list(parse(raw))
        assert len(data) == count
        for entry in data:
            assert 'TPA:' not in entry.description
            assert entry.description
            for reference in entry.references:
                assert reference.title != ';'


def test_creates_simple_entry():
    with open('data/ena/ncr/wgs/aa/wgs_aacd01_fun.ncr', 'rb') as raw:
        simple_data = list(parse(raw))

    assert attr.asdict(simple_data[0]) == attr.asdict(Entry(
        primary_id='',
        accession='AACD01000002.1:101667..101773:tRNA',
        ncbi_tax_id=227321,
        database='ENA',
        sequence='GCCCGGATGGTCTAGTGGTATGATTCTCCCTTCGGGGGCAGCGCCCGGTACATAATAACATGTATCAGAAATGGGAGAGGTCCCGCGTTCGAATCGCGGTTCGGGCC',
        exons=[Exon(
            chromosome='VIII',
            complement=False,
            primary_end=101773,
            primary_start=101667,
        )],
        rna_type='tRNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AACD01000002.1:101667..101773:tRNA',
        seq_version='1',
        xref_data={'ena_refs': {'BIOSAMPLE': ('SAMN02953587', None)}},
        chromosome='VIII',
        species='Aspergillus nidulans FGSC A4',
        lineage=(
            'Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; '
            'Eurotiomycetes; Eurotiomycetidae; Eurotiales; Aspergillaceae; '
            'Aspergillus; Aspergillus nidulans FGSC A4'
        ),
        product='tRNA-Pro',
        parent_accession='AACD01000002',
        project='PRJNA130',
        keywords='WGS',
        division=None,
        # division='FUN',
        description='Aspergillus nidulans FGSC A4 tRNA-Pro',
        mol_type='genomic DNA',
        is_composite='N',
        references=[
            Reference(
                accession='AACD01000002.1:101667..101773:tRNA',
                authors=(
                    "Galagan J.E., Calvo S.E., Cuomo C., Ma L.J., Wortman J.R."
                    ", Batzoglou S., Lee S.I., Basturkmen M., Spevak C.C., "
                    "Clutterbuck J., Kapitonov V., Jurka J., Scazzocchio C., "
                    "Farman M., Butler J., Purcell S., Harris S., Braus G.H., "
                    "Draht O., Busch S., D'Enfert C., Bouchier C., Goldman "
                    "G.H., Bell-Pedersen D., Griffiths-Jones S., Doonan J.H., "
                    "Yu J., Vienken K., Pain A., Freitag M., Selker E.U., "
                    "Archer D.B., Penalva M.A., Oakley B.R., Momany M., "
                    "Tanaka T., Kumagai T., Asai K., Machida M., Nierman W.C.,"
                    " Denning D.W., Caddick M., Hynes M., Paoletti M., "
                    "Fischer R., Miller B., Dyer P., Sachs M.S., Osmani S.A., "
                    "Birren B.W."
                ),
                location='Nature 438(7071):1105-1115(2005).',
                title=(
                    'Sequencing of Aspergillus nidulans and comparative '
                    'analysis with A. fumigatus and A. oryzae'
                ),
                pmid=16372000,
                doi=None,
            ),
            Reference(
                accession='AACD01000002.1:101667..101773:tRNA',
                location=(
                    'Submitted (04-APR-2003) to the INSDC. Whitehead '
                    'Institute/MIT Center for Genome Research, 320 Charles '
                    "Street, Cambridge, MA 02142, USA"
                ),
                title=None,
                authors=(
                    'Birren B., Nusbaum C., Abouelleil A., Allen N., '
                    'Anderson S., Arachchi H.M., Barna N., Bastien V., '
                    'Bloom T., Boguslavkiy L., Boukhgalter B., Butler J., '
                    'Calvo S.E., Camarata J., Chang J., Choepel Y., '
                    'Collymore A., Cook A., Cooke P., Corum B., DeArellano K.,'
                    ' Diaz J.S., Dodge S., Dooley K., Dorris L., Elkins T., '
                    'Engels R., Erickson J., Faro S., Ferreira P., '
                    'FitzGerald M., Gage D., Galagan J., Gardyna S., Gnerre S.'
                    ', Graham L., Grand-Pierre N., Hafez N., Hagopian D., '
                    'Hagos B., Hall J., Horton L., Hulme W., Iliev I., '
                    'Jaffe D., Johnson R., Jones C., Kamal M., Kamat A., '
                    'Karatas A., Kells C., Landers T., Levine R., '
                    'Lindblad-Toh K., Liu G., Lui A., Ma L.-J., Mabbitt R., '
                    'MacLean C., Macdonald P., Major J., Manning J., '
                    'Matthews C., Mauceli E., McCarthy M., Meldrim J., '
                    'Meneus L., Mihova T., Mlenga V., Murphy T., Naylor J., '
                    "Nguyen C., Nicol R., Nielsen C.B., Norbu C., O'Connor T.,"
                    " O'Donnell P., O'Neil D., Oliver J., Peterson K., "
                    'Phunkhang P., Pierre N., Purcell S., Rachupka A., '
                    'Ramasamy U., Raymond C., Retta R., Rise C., Rogov P., '
                    'Roman J., Schauer S., Schupback R., Seaman S., Severy P.,'
                    ' Smirnov S., Smith C., Spencer B., Stange-Thomann N., '
                    'Stojanovic N., Stubbs M., Talamas J., Tesfaye S., '
                    'Theodore J., Topham K., Travers M., Vassiliev H., '
                    'Venkataraman V.S., Viel R., Vo A., Wang S., Wilson B., '
                    'Wu X., Wyman D., Young G., Zainoun J., Zembek L., '
                    'Zimmer A., Zody M., Lander E.'
                ),
                pmid=None,
                doi=None,
            ),
            Reference(
                accession='AACD01000002.1:101667..101773:tRNA',
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
                    'Submitted (07-JAN-2004) to the INSDC. Whitehead '
                    'Institute/MIT Center for Genome Research, '
                    '320 Charles Street, Cambridge, MA 02142, USA'
                ),
                title=None,
                pmid=None,
                doi=None,
            )
        ]
    ))


def test_can_find_correct_ncRNA_type():
    with open('data/ena/ncr/wgs/aa/wgs_abxv02_pro.ncr', 'rb') as raw:
        ncrna_data = list(parse(raw))

    assert attr.asdict(ncrna_data[0]) == attr.asdict(Entry(
        primary_id='',
        accession='ABXV02000002.1:33498..33573:ncRNA',
        ncbi_tax_id=500637,
        database='ENA',
        sequence='ACTGCTTTTCTTTGATGTCCCCATATTGAGGAGCCCGATAGCCATTTGATTACTTCATGCTATCGGGTTTTTTATT',
        exons=[Exon(
            chromosome='',
            complement=False,
            primary_end=33573,
            primary_start=33498
        )],
        rna_type='other',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:ABXV02000002.1:33498..33573:ncRNA',
        seq_version='1',
        note_data={
            'text': ['Rfam score 66']
        },
        xref_data={'ena_refs': {'BIOSAMPLE': ('SAMN00008855', None)}},
        species="Providencia rustigianii DSM 4541",
        lineage=(
            "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; "
            "Morganellaceae; Providencia; Providencia rustigianii DSM 4541"
        ),
        product="RybB RNA",
        parent_accession='ABXV02000002',
        keywords='WGS',
        division=None,
        # division='PRO',
        description='Providencia rustigianii DSM 4541 RybB RNA',
        mol_type="genomic DNA",
        is_composite='N',
        inference="nucleotide motif:Rfam:RF00110",
        locus_tag="PROVRUST_04548",
        project='PRJNA28651',
        references=[
            Reference(
                accession='ABXV02000002.1:33498..33573:ncRNA',
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
                accession='ABXV02000002.1:33498..33573:ncRNA',
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
                accession='ABXV02000002.1:33498..33573:ncRNA',
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
    ))


def test_can_parse_all_example_entries():
    with open('data/test_example_entries.ncr', 'rb') as raw:
        examples = list(parse(raw))

    assert attr.asdict(examples[0]) == attr.asdict(Entry(
        primary_id='',
        accession='AB330787.1:1..34:misc_RNA',
        ncbi_tax_id=9606,
        database='ENA',
        sequence="ATTGGGGAGTGAGAAGGAGAGAACGCGGTCTGAA",
        exons=[Exon(
            chromosome='',
            complement=False,
            primary_start=1,
            primary_end=34,
        )],
        rna_type='misc_RNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330787.1:1..34:misc_RNA',
        seq_version='1',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division=None,
        # division='HUM',
        description='Homo sapiens (human) miscellaneous RNA',
        parent_accession='AB330787',
        note_data={
            'text': [
                "RNA sequence binds to tRNaseZL",
                "similar to U3 snRNA fragment"
            ]
        },
        mol_type='other RNA',
        is_composite='N',
    ))

    assert attr.asdict(examples[1]) == attr.asdict(Entry(
        primary_id='',
        accession='AB330786.1:1..27:misc_RNA',
        database='ENA',
        ncbi_tax_id=9606,
        sequence='ATTGCAGTACCTCCAGGAACGGTGCAC',
        exons=[Exon(
            chromosome='',
            complement=False,
            primary_start=1,
            primary_end=27,
        )],
        rna_type='misc_RNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330786.1:1..27:misc_RNA',
        seq_version='1',
        xref_data={'ena_refs': {'RFAM': ('RF00005', 'tRNA')}},
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division=None,
        # division='HUM',
        description='Homo sapiens (human) miscellaneous RNA',
        parent_accession='AB330786',
        note_data={
            'text': [
                'RNA sequence binds to tRNaseZL',
                'similar to U2 snRNA fragment'
            ]
        },
        mol_type='other RNA',
        is_composite='N',
    ))

    assert attr.asdict(examples[2]) == attr.asdict(Entry(
        primary_id='',
        accession='AB330785.1:1..34:misc_RNA',
        database='ENA',
        ncbi_tax_id=9606,
        sequence='CGCGACCTCAGATCAGACGTGGCGACCCGCTGAA',
        exons=[Exon(
            chromosome='',
            complement=False,
            primary_start=1,
            primary_end=34,
        )],
        rna_type='misc_RNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330785.1:1..34:misc_RNA',
        seq_version='1',
        xref_data={'ena_refs': {'SRPDB': ('Xeno.laev._DC015625', None)}},
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division=None,
        # division='HUM',
        description='Homo sapiens (human) miscellaneous RNA',
        note_data={
            'text': [
                "RNA sequence binds to tRNaseZL",
                "similar to 28S rRNA fragment",
            ]
        },
        parent_accession='AB330785',
        mol_type='other RNA',
        is_composite='N',
        keywords='RNAcentral; Third Party Annotation; TPA; TPA:specialist_db',
    ))

    assert attr.asdict(examples[3]) == attr.asdict(Entry(
        primary_id='',
        accession='HAAO01001079.1:1..21:ncRNA',
        database='ENA',
        ncbi_tax_id=9606,
        sequence='ACCACTGCACTCCAGCCTGAG',
        exons=[Exon(
            chromosome='',
            complement=False,
            primary_start=1,
            primary_end=21,
        )],
        rna_type='miRNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:HAAO01001079.1:1..21:ncRNA',
        seq_version='1',
        xref_data={'ena_refs': {'MIRBASE': ('MIMAT0022742', 'hsa-miR-1273g-3p')}},
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division=None,
        # division='HUM',
        project='PRJEB4451',
        description='Homo sapiens (human) microRNA hsa-miR-1273g-3p',
        product="microRNA hsa-miR-1273g-3p",
        experiment="EXISTENCE:RNA-seq ECO0000205",
        mol_type='transcribed RNA',
        gene="hsa-miR-1273g-3p",
        keywords=(
            'RNAcentral; TPA; TPA:specialist_db; Transcriptome Shotgun '
            'Assembly; TSA'
        ),
        is_composite='N',
        parent_accession='HAAO01001079',
        references=[
            Reference(
                accession='HAAO01001079.1:1..21:ncRNA',
                authors='',
                location='Submitted (20-AUG-2013) to the INSDC.',
                title=None,
                pmid=None,
                doi=None,
            ),
            Reference(
                accession='HAAO01001079.1:1..21:ncRNA',
                authors='Kozomara A., Griffiths-Jones S.',
                location='Nucleic Acids Res. 39(Database issue):D152-D157(2011).',
                title='miRBase: integrating microRNA annotation and deep-sequencing data',
                pmid=21037258,
                doi=None,
            ),
        ]
    ))

    assert attr.asdict(examples[4]) == attr.asdict(Entry(
        primary_id='',
        accession='HG519048.1:1..359:tmRNA',
        database='ENA',
        ncbi_tax_id=32644,
        sequence="GGGAGCGACTTGGCTTCGACAGGAGTAAGTCTGCTTAGATGGCATGTCGCTTTGGGCAAAGCGTAAAAAGCCCAAATAAAATTAAACGCAAACAACGTTAAATTCGCTCCTGCTTACGCTAAAGCTGCGTAAGTTCAGTTGAGCCTGAAATTTAAGTCATACTATCTAGCTTAATTTTCGGTCATTTTTGATAGTGTAGCCTTGCGTTTGACAAGCGTTGAGGTGAAATAAAGTCTTAGCCTTGCTTTTGAGTTTTGGAAGATGAGCGAAGTAGGGTGAAGTAGTCATCTTTGCTAAGCATGTAGAGGTCTTTGTGGGATTATTTTTGGACAGGGGTTCGATTCCCCTCGCTTCCACCA",
        exons=[Exon(
            chromosome='',
            complement=False,
            primary_start=1,
            primary_end=359,
        )],
        rna_type='tmRNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:HG519048.1:1..359:tmRNA',
        seq_version='1',
        note_data={
            'text': ['Tag:(A)ANNVKFAPAYAKAA*'],
        },
        xref_data={'ena_refs': {'TMRNA-WEBSITE': ('Campy_jejun_700819', None)}},
        lineage='unclassified sequences; unidentified',
        species='unidentified',
        division=None,
        # division='UNC',
        keywords='RNAcentral; TPA; TPA:specialist_db',
        description='unidentified transfer-messenger mRNA Campy_jejun_700819',
        product="transfer-messenger mRNA Campy_jejun_700819",
        gene="tmRNA Campy_jejun_700819",
        project='PRJEB4570',
        mol_type='genomic DNA',
        parent_accession='HG519048',
        is_composite='N',
        references=[
            Reference(
                accession='HG519048.1:1..359:tmRNA',
                authors='',
                location='Submitted (13-SEP-2013) to the INSDC.',
                title=None,
                pmid=None,
                doi=None,
            ),

            Reference(
                accession='HG519048.1:1..359:tmRNA',
                authors='Gueneau de Novoa P., Williams K.P.',
                location='Nucleic Acids Res. 32(Database issue):D104-D108(2004).',
                title=(
                    "The tmRNA website: reductive evolution of tmRNA in "
                    "plastids and other endosymbionts"
                ),
                pmid=14681369,
                doi=None,
            ),
        ]
    ))

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
    #     rna_type='snoRNA',
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
    #     rna_type='misc_RNA'
    #     mol_type='genomic DNA',
    #     locus_tag="SPNCRNA.1210",
    #     chromosome='III',
    #     product='antisense RNA (predicted)'
    # )),


def test_can_handle_file_with_invalid_fields():
    with open('data/test_invalid_fields.ncr', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'KM079256.1:1..1300:rRNA'
    assert data.mol_type == 'genomic DNA'
    assert data.species == 'Candidatus Stammerula sp. of Trupanea "pohakuloa"'
    assert data.ncbi_tax_id == 1630665
    assert data.product == "16S ribosomal RNA"


def test_can_parse_out_anticodon_from_gene():
    with open('data/ena/trna-with-anticodon.embl', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'HF536610.1:1288..1686:tRNA'
    assert data.ncbi_tax_id == 188169
    assert data.product == "transfer RNA Leu"
    assert data.gene == "tRNA-Leu (UAA)"
    assert data.mol_type == 'genomic DNA'
    assert data.anticodon == 'UAA'
    assert data.organelle == 'plastid:chloroplast'


def test_can_parse_anticodon_from_gtrnadb_stype_gene():
    with open('data/ena/trna-with-anticodon.embl', 'rb') as raw:
        data = list(parse(raw))[1]

    assert data.accession == 'LK008175.1:1..73:tRNA'
    assert data.ncbi_tax_id == 411154
    assert data.mol_type == "genomic DNA"
    assert data.anticodon == 'GCC'
    assert data.product == "transfer RNA-Gly (GCC)"
    assert data.gene == "tRNA-Gly-GCC-1-1"
    assert data.standard_name == "tRNA-Gly (GCC)"
    assert data.inference == "ab initio prediction:tRNAscan-SE:1.21"
    assert data.note_data == {
        'ontology': [
            "ECO:0000202",
            "GO:0030533",
            "SO:0000253",
        ],
        "text": [
            "Covariance Model: Bacteria; CM Score: 87.61",
            "Legacy ID: chr.trna3-GlyGCC"
        ]
    }


def test_can_parse_anticodon_from_note():
    with open('data/ena/anticodon-in-note.embl', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'CP000102.1:337323..337406:tRNA'
    assert data.note_data == {'text': ['codon recognized: CUA']}
    assert data.anticodon == 'UAG'
    assert data.operon == 'trnA'
    assert data.product == 'tRNA-Leu'
    assert data.gene is None
    assert data.locus_tag == 'Msp_0274'


def test_can_parse_pseudogene():
    with open('data/ena/pseudogene.embl', 'rb') as raw:
        data = next(parse(raw))
    assert data.accession == 'NIDN01000248.1:29758..29889:tRNA'
    assert data.pseudogene == 'unprocessed'


def test_can_parse_function():
    with open('data/ena/function.embl', 'rb') as raw:
        data = list(parse(raw))

    assert data[0].accession == 'EU410654.1:1..92:ncRNA'
    assert data[0].function == 'guide for 26S rRNA methylation at U1043'
    assert data[0].references[0].pmid == 18493037
    # TODO: Improve biopython so it can parse DOI's from references
    # assert data[0].references[0].doi == '10.1534/genetics.107.086025'

    assert data[1].accession == 'AB046489.1:221..306:tRNA'
    assert data[1].function == 'tRNA-Pro'
    assert data[1].organelle == 'mitochondrion'

    assert data[2].accession == 'CP003783.1:1548698..1548818:ncRNA'
    assert data[2].function == '1.8: Sporulation'
    assert data[2].gene == "csfG"
    assert data[2].product == "sporulation-specific regulatory RNA"
    assert data[2].locus_tag == "B657_miscRNA23"
    assert data[2].rna_type == 'other'


def test_can_parse_old_locus_tag():
    with open('data/ena/old_locus_tag.embl', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'AL591985.1:1315182..1315258:tRNA'
    assert data.parent_accession == 'AL591985'
    assert data.project == 'PRJNA19'
    assert data.old_locus_tag == 'SMb21712'
    assert data.locus_tag == "SM_b21712"
    assert data.gene == "tRNA-ARG_CCG"
    assert data.anticodon == "CCG"


def test_can_parse_operons():
    with open('data/ena/operons.embl', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'CP000102.1:1163547..1163662:rRNA'
    assert data.operon == 'rrnD'


def test_can_parse_gene_synonyms():
    with open('data/ena/gene_synonym.embl', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'CP000948.1:2011661..2011909:misc_RNA'
    assert data.gene_synonyms == ['IS091', 'sraC', 'tpke79']
    assert data.locus_tag == 'ECDH10B_1978'
    assert data.mol_type == 'genomic DNA'
    assert data.gene == 'ryeA'
    assert data.product == 'small RNA'
    assert data.project == 'PRJNA20079'


@pytest.mark.parametrize('pmid', [
    12244299, 6196367, 6181418, 6802847, 3403542, 6084597, 10924331, 7528809,
    7529207, 1704372, 20610725, 18617187, 17881443, 17164479, 8389475,
    10834842, 10684931, 15611297, 20668672, 911771, 6209580])
def test_can_extract_references_from_experiment(pmid):
    with open('data/ena/experiment-references.embl', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'HG975378.1:1..299:ncRNA'
    assert data.rna_type == 'SRP_RNA'
    assert data.experiment == 'EXISTENCE: lncRNAdb literature review [PMID: 12244299,6196367,6181418,6802847,3403542,6084597,10924331, 7528809,7529207,1704372,20610725,18617187,17881443, 17164479,8389475,10834842,10684931,15611297,20668672, 911771,6209580]'
    assert data.gene == 'RN7SL1'
    assert data.mol_type == 'transcribed RNA'
    assert data.product == 'Small nucleolar RNA 7SL'
    assert data.note_data == {
        'ontology': ['ECO:0000305', 'GO:0006617', 'GO:0048501', 'SO:0000590'],
        'text': ['biotype:SRP_RNA'],
    }
    assert pmid in [ref.pmid for ref in data.references]


@pytest.mark.parametrize('pmid', [
    1379177, 8288542, 9098041, 9826762, 12165569, 12547201, 1317842, 15831787
])
def test_can_extract_references_from_note_or_experiment(pmid):
    with open('data/ena/note-references.embl', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'AL009126.3:3856226..3856447:misc_RNA'
    assert data.note_data == {
        'text': ['Evidence 1a: Function experimentally demonstrated in the studied strain; PubMedId: 1379177, 8288542, 9098041, 9826762, 12165569, 12547201, 1317842, 15831787; Product type n: RNA']
    }
    assert data.product == 'T-box'
    assert data.locus_tag == 'BSU_misc_RNA_57'
    assert data.experiment == 'publication(s) with functional evidences, PMID:1379177, 8288542, 9098041, 9826762, 12165569, 12547201, 1317842, 15831787'
    assert data.function == '16.3: Control'
    assert data.gene == 'tboTB'
    assert pmid in [ref.pmid for ref in data.references]


def test_can_handle_unclosed_parens():
    with open('data/test_feature_unclosed_parenthesis.ncr', 'rb') as raw:
        data = next(parse(raw))

    assert data.accession == 'HE860504.1:1..14644:tRNA'
    assert data.mol_type == 'genomic DNA'
    assert data.ncbi_tax_id == 1200666
    assert data.gene == 'tRNA-Ser'
    assert data.product == "transfer RNA Serine"
    assert data.species == "Metacrangonyx sp. 3 ssp. 1 MDMBR-2012"
    assert data.organelle == "mitochondrion"
    assert data.project is None
    assert data.description == 'Metacrangonyx sp. 3 ssp. 1 MDMBR-2012 transfer RNA Serine'
    assert data.anticodon == '(pos:HE860504.1: complement(14629..14631),aa:Ser)'


def test_can_extract_xrefs():
    with open('data/test_example_entries.ncr', 'rb') as raw:
        data = list(parse(raw))

    assert data[5].accession == 'HG975377.1:1..332:ncRNA'
    assert data[5].xref_data == {'lncrnadb': ['190', '7SK']}
    assert data[7].accession == 'LM608264.1:7..26:ncRNA'
    assert data[7].xref_data == {'mirbase': ['MI0000182']}
