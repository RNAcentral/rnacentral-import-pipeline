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

import pytest

from rnacentral_pipeline.rnacentral.precompute.description import \
    description_of

from .helpers import load_data


@pytest.mark.parametrize('rna_id,rna_type,name', [  # pylint: disable=no-member
    ('URS000001E7BA_559292', 'tRNA', 'Saccharomyces cerevisiae S288c tRNA-Gln (tQ(UUG)C, tQ(UUG)D1-3, tQ(UUG)E1, tQ(UUG)H, tQ(UUG)L)'),
    ('URS00000AEE53_380749', 'tmRNA', 'Hydrogenobaculum sp. Y04AAS1 tmRNA'),
    ('URS00000F9D45_9606', 'rRNA', 'Homo sapiens RNA, 5S ribosomal 1 (RNA5S1-8, RNA5S10-17)'),
    ('URS000018EB2E_3702', 'lncRNA', 'Arabidopsis thaliana (thale cress) Long non-coding antisense RNA COOLAIR'),
    ('URS000019E0CD_9606', 'lncRNA', 'Homo sapiens HELLP associated long non-coding RNA (HELLPAR)'),
    ('URS00001DEEBE_562', 'tRNA', 'Escherichia coli tRNA-Pro (CGG) (tRNA-Pro-CGG-1-1)'),
    ('URS00002F21DA_7227', 'precursor_RNA', 'Drosophila melanogaster bantam stem-loop (dme-bantam)'),
    ('URS000034C5CB_7227', 'SRP_RNA', 'Drosophila melanogaster (fruit fly) Signal recognition particle 7SL RNA CR32864 (Dmel_CR32864, Dmel_CR42652)'),
    ('URS000037602E_9606', 'tmRNA', 'Homo sapiens (human) transfer-messenger RNA Esche_coli_K12'),
    ('URS00003AC4AA_3702', 'siRNA', 'Arabidopsis thaliana (thale-cress) tAS3/TASIR-ARF (TRANS-ACTING SIRNA3); other RNA'),
    ('URS00003BECAC_9606', 'lncRNA', 'Homo sapiens (human) long intergenic non-protein coding RNA 1729 (LINC01729)'),
    ('URS00003CE153_9606', 'lncRNA', 'Homo sapiens STARD4 antisense RNA 1 (STARD4-AS1)'),
    ('URS00003EBD9A_9913', 'telomerase_RNA', 'Bos taurus telomerase RNA component (TERC)'),
    ('URS0000466DE6_6239', 'miRNA', 'Caenorhabditis elegans cel-miR-229-5p'),
    ('URS000048B30C_3702', 'tRNA', 'Arabidopsis thaliana (thale cress) partial tRNA-Leu'),
    ('URS00004E52D3_10090', 'lncRNA', 'Mus musculus (house mouse) predicted gene 12238 (Gm12238)'),
    ('URS00004E9E38_7227', 'miRNA', 'Drosophila melanogaster (fruit fly) dme-bantam-3p'),
    ('URS00004FB44B_6239', 'rRNA', 'Caenorhabditis elegans 26s rRNA'),
    pytest.param('URS000051DCEC_10090', 'snoRNA', 'Mus musculus small nucleolar RNA, C/D box 17 (Snord17), small nucleolar RNA',
                 marks=pytest.mark.xfail),
    ('URS00005511ED_6239', 'lncRNA', 'Caenorhabditis elegans long non-coding RNA linc-125'),
    ('URS000055786A_7227', 'miRNA', 'Drosophila melanogaster (fruit fly) dme-bantam-5p'),
    ('URS0000563A36_7227', 'snoRNA', 'Drosophila melanogaster (fruit fly) snoRNA:Tudor-SN-a (Dmel_CR43585)' ),
    ('URS0000569A4A_9606', 'snoRNA', 'Homo sapiens (human) small Cajal body-specific RNA 10 (SCARNA10)'),
    ('URS00005F4CAF_3702', 'tRNA', 'Arabidopsis thaliana (thale cress) tRNA-Met(CAT)'),
    pytest.param('URS000060B496_10090', 'snoRNA', 'Mus musculus small nucleolar RNA, H/ACA box 3 (Snora3), small nucleolar RNA',
                 marks=pytest.mark.xfail),
    ('URS000061F377_559292', 'rRNA', 'Saccharomyces cerevisiae S288c (RDN25-1, RDN25-2)'),
    ('URS00006550DA_10090', 'snoRNA', 'Mus musculus (house mouse) small Cajal body-specific RNA 1 (Scarna13)'),
    ('URS0000661037_7955', 'tRNA', 'Danio rerio tRNA'),
    ('URS000069D7FA_6239', 'tRNA', 'Caenorhabditis elegans tRNA-His'),
    ('URS00006B3271_10090', 'snoRNA', 'Mus musculus (house mouse) small Cajal body-specific RNA 2 (Scarna2)'),
    ('URS00006CE02F_9606', 'snoRNA', 'Homo sapiens (human) snoRNA U8 small nucleolar RNA (RF00096)'),
    ('URS00006D80BC_9913', 'precursor_RNA', 'Bos taurus (cattle) microRNA bta-mir-497 precursor'),
    ('URS00006DC8B9_6239', 'tRNA', 'Caenorhabditis elegans tRNA-Undet'),
    ('URS00007150F8_9913', 'precursor_RNA', 'Bos taurus (cattle) microRNA bta-mir-431 precursor'),
    ('URS0000759BEC_9606', 'lncRNA', 'Homo sapiens DiGeorge syndrome critical region gene 9 (DGCR9)'),
    ('URS000075A3E2_9606', 'precursor_RNA', 'Homo sapiens (human) microRNA hsa-mir-6859 precursor (hsa-mir-6859 1 to 4)'),
    ('URS000075A546_9606', 'precursor_RNA', 'Homo sapiens (human) microRNA hsa-mir-3648 precursor (hsa-mir-3648-1, hsa-mir-3648-2)'),
    ('URS000075B196_7955', 'precursor_RNA', 'Danio rerio (zebrafish) microRNA dre-mir-430c precursor (dre-mir-430c 1 to 18)'),
    ('URS000075C808_9606', 'lncRNA', 'Homo sapiens HOX transcript antisense RNA (HOTAIR)'),
    ('URS000075CC93_9606', 'precursor_RNA', 'Homo sapiens (human) microRNA hsa-mir-1302 precursor (hsa-mir-1302-2, hsa-mir-1302 9 to 11)'),
    ('URS000075CF25_9913', 'precursor_RNA', 'Bos taurus (cattle) microRNA bta-mir-10a precursor'),
    ('URS0000795103_6239', 'precursor_RNA', 'Caenorhabditis elegans microRNA cel-mir-8204 precursor'),
    pytest.param('URS0000808D70_1478174', 'tmRNA', 'Neochlamydia sp. TUME1 tmRNA',
                 marks=pytest.mark.xfail(reason="Inactive")),
    ('URS00008E3A1B_10090', 'lncRNA', 'Mus musculus predicted gene 11532 (Gm11532)'),
    ('URS00009E8F92_885695', 'rRNA', 'Sinumelon nullarboricum partial 16S ribosomal RNA'),
    pytest.param('URS0000A767C0_3702', 'lncRNA', 'Arabidopsis thaliana (thale-cress) potential natural antisense gene, locus overlaps with AT1G44120',
                 marks=pytest.mark.xfail),
    ('URS0000A86584_10090', 'ncRNA', 'Mus musculus (house mouse) predicted gene 29254 (Gm29254)'),
    ('URS0000ABD87F_9606', 'rRNA', 'Homo sapiens RNA, 45S pre-ribosomal N1 (RNA45SN1)'),
    ('URS0000D50284_7240', 'precursor_RNA', 'Drosophila simulans microRNA dsi-mir-988 precursor'),
    ('URS0000D5657E_7240', 'precursor_RNA', 'Drosophila simulans microRNA dsi-mir-4966 precursor (dsi-mir-4966-3)'),
    ('URS00003CF845_9606', 'miRNA', 'Homo sapiens (human) hsa-miR-1273d'),
    pytest.param('URS0000BD9FEB_10089', 'ribozyme', 'Mus caroli (Ryukyu mouse) nuclear RNase P',
                 marks=pytest.mark.xfail),

    # TODO: Implement the logic required below
    # ('URS000041AF00_274', 'tRNA', 'TRNA-FMET from Thermus thermophilus (29 structures)'),
    # ('URS0000A7633F_9606', 'snRNA', "U4 snRNA from Homo sapiens (PDB 5H1K, chain C,D"),
    # ('URS0000D2348C_9606', 'misc_RNA', 'Homo sapiens (human) let-7 microRNA precursor (2 structures)'),
    # ('URS000080DD8D_32630', 'misc_RNA', "25-nt RNA from synthetic construct (PDB 4V8X, chain AX)"),
    # ('URS000080DE2B_274', 'rRNA', "5S rRNA from Thermus thermophilus (34 structures)"),
    # ('URS00000C777C_1231464', "Alkalibacterium sp. 3.5P*23 small subunit ribosomal RNA (16S)"),
    # ('URS000080DD59_32630', "5'-R(*UP*GP*(CBV)P*(CBV)P*AP*GP*UP*UP*CP*GP*CP*UP*GP*GP*C)-3' from (PDB 1QBP, chain E)"),
    # ('URS000080DE07_32630', "RNA (5'-R(P*CP*GP*AP*UP*CP*GP*GP*GP*UP*GP*UP*C)-3') from (PDB 1RY1, chain Q)"),
    # ('URS000080DE5B_32630', "RNA (5'-R(P*AP*UP*CP*GP*CP*GP*CP*CP*UP*GP*UP*G)-3') from (PDB 1RY1, chain R)"),
    # ('URS000080E209_32630', "RNA (5'-R(*CP*CP*GP*CP*CP*GP*CP*GP*CP*CP*AP*(5BU)P*GP*CP*CP*UP*GP*UP*GP*GP*CP*GP... from (PDB 3MEI, chain B)"),
    # ('URS000080E22B_308052', "tRNA (5'-D(*AP*UP*CP*CP*CP*CP*GP*UP*GP*UP*CP*CP*UP*UP*GP*GP*UP*UP*CP*G)-3') from Mitsuaria sp. 67 (PDB 4WT8, chain D2)"),
    # ('URS000080E230_32630', "5'-D(*CP*AP*GP*CP*TP*AP*CP*TP*TP*GP*AP*GP*CP*T)-3' from (PDB 3H3V, chain P)"),
    # ('URS0000A77852_32630', "RNA (5'-R(*(LCC)P*(LCC)P*(LCA)P*(LCG)P*AP*CP*UP*UP*AP*AP*GP*UP*CP*U)-3') from synthetic construct (PDB 5V0J, chain B)"),
    # ('URS0000BC4693_9606', "RNA (5'-R(*GP*(CBV)P*CP*GP*GP*(6MZ)P*UP*GP*GP*C)-3') from Homo sapiens (PDB 5LR4, chain D)"),
    # ('URS000080E10B_32630', "5'-R(*GP*CP*CP*GP*AP*AP*GP*CP*CP*(P5P)-3' from (PDB 1XV0, chain B)"),
    # ('URS000080E135_274', "messenger RNA (5'-R(*AP*AP*UP*GP*UP*AP*G)-3') from Thermus thermophilus (PDB 4V9N, chain CV)"),
    # ('URS000041AF00_274', 'RNA (77-MER) from Thermus thermophilus (PDB 4V7J, chain Bw)'),
])
def test_computes_correct_species_specific_descriptions(rna_id, rna_type, name):
    data = load_data(rna_id)
    assert name == description_of(rna_type, data)


# @pytest.mark.parametrize('upi,rna_type,name', [
#     ('URS0000181AEC', 'tRNA', 'tRNA from 3413 species'),
# ])
# def test_computes_correct_generic_names(upi, rna_type, name):
#     data = fetch.entry(db(), upi)
#     assert description_of(rna_type, data) == name
