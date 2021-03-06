# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.europepmc import fetch

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.data import IdReference


def lookup(ref_id):
    return attr.asdict(fetch.lookup(IdReference.build(ref_id)))


@pytest.mark.parametrize('raw_id', [
    28815543,
    '28815543',
    'PMC5890441',
    'doi:10.1007/978-981-10-5203-3_9',
])
def test_can_fetch_publication(raw_id):
    assert fetch.summary(IdReference.build(raw_id)) == {
        "id": "28815543",
        "source": "MED",
        "pmid": "28815543",
        "doi": "10.1007/978-981-10-5203-3_9",
        "title": "Understanding the Role of lncRNAs in Nervous System Development.",
        "authorString": "Clark BS, Blackshaw S.",
        "journalTitle": "Adv Exp Med Biol",
        "journalVolume": "1008",
        "pubYear": "2017",
        "journalIssn": "0065-2598; 2214-8019; ",
        "pageInfo": "253-282",
        "pubType": "research-article; review; journal article; research support, n.i.h., extramural",
        "isOpenAccess": "N",
        "inEPMC": "Y",
        "inPMC": "N",
        "hasPDF": "Y",
        "hasBook": "N",
        "citedByCount": 4,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2017-01-01",
        'hasSuppl': 'N',
        'pmcid': 'PMC5890441',
        'firstIndexDate': '2017-08-18',
    }


def test_complains_given_bad_pmid():
    with pytest.raises(Exception):
        fetch.summary(IdReference.build(-1))


@pytest.mark.parametrize('raw_id', [
    27858507,
    '27858507',
    'doi:10.1080/15476286.2016.1251002',
    'PMCID:PMC5785218',
])
def test_can_build_reference(raw_id):
    assert lookup(raw_id) == attr.asdict(Reference(
        authors=u"Lünse CE, Weinberg Z, Weinberg Z, Breaker RR.",
        location='RNA Biol 14(11):1499-1507 (2017)',
        title=(
            'Numerous small hammerhead ribozyme variants associated with '
            'Penelope-like retrotransposons cleave RNA as dimers'
        ),
        pmid=27858507,
        doi='10.1080/15476286.2016.1251002',
        pmcid='PMC5785218',
    ))


def test_can_deal_with_unicode():
    data = fetch.lookup(IdReference.build(27334534))
    assert attr.asdict(data) == {
        'authors': 'Ahn S, Jin TE, Chang DH, Rhee MS, Kim HJ, Lee SJ, Park DS, Kim BC.',
        'location': 'Int J Syst Evol Microbiol 66(9):3656-3661 (2016)',
        'title': 'Agathobaculum butyriciproducens gen. nov. \xa0sp. nov., a strict anaerobic, butyrate-producing gut bacterium isolated from human faeces and reclassification of Eubacterium desmolans as Agathobaculum desmolans comb. nov',
        'pmid': 27334534,
        'doi': '10.1099/ijsem.0.001195',
        'pmcid': None,
    }


def test_builds_correction_location():
    assert lookup(26184978) == attr.asdict(Reference(
        authors='Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.',
        location='Sci Rep 5:12276 (2015)',
        title=(
            'MiR-135b-5p and MiR-499a-3p Promote Cell '
            'Proliferation and Migration in Atherosclerosis by Directly '
            'Targeting MEF2C'
        ),
        pmid=26184978,
        doi='10.1038/srep12276',
        pmcid='PMC4505325',
    ))


def test_can_handle_missing_volume():
    assert lookup(27389411) == attr.asdict(Reference(
        authors='Carr G, Barrese V, Stott JB, Povstyan OV, Jepps TA, Figueiredo HB, Zheng D, Jamshidi Y, Greenwood IA.',
        location='Cardiovasc Res 112(2):581-589 (2016)',
        title='MicroRNA-153 targeting of KCNQ4 contributes to vascular dysfunction in hypertension',
        pmid=27389411,
        doi='10.1093/cvr/cvw177',
        pmcid='PMC5079273',
    ))


def test_it_can_find_if_duplicate_ext_ids():
    assert lookup(375006) == attr.asdict(Reference(
        authors='Macino G, Tzagoloff A.',
        location='Mol Gen Genet 169(2):183-188 (1979)',
        title='Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae',
        pmid=375006,
        doi='10.1007/bf00271669',
        pmcid=None,
    ))


def test_can_lookup_by_doi():
    assert lookup('doi:10.1007/bf00271669') == attr.asdict(Reference(
        authors='Macino G, Tzagoloff A.',
        location='Mol Gen Genet 169(2):183-188 (1979)',
        title='Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae',
        pmid=375006,
        doi='10.1007/bf00271669',
        pmcid=None,
    ))


@pytest.mark.parametrize('ref_id', [
    26184978,
    'pmid:26184978',
    'PMID:26184978',
    'doi:10.1038/srep12276',
    'DOI:10.1038/srep12276',
])
def test_can_handle_several_reference_formats(ref_id):
    assert lookup(ref_id) == attr.asdict(Reference(
        authors='Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.',
        location='Sci Rep 5:12276 (2015)',
        title=(
            'MiR-135b-5p and MiR-499a-3p Promote Cell '
            'Proliferation and Migration in Atherosclerosis by Directly '
            'Targeting MEF2C'
        ),
        pmid=26184978,
        doi='10.1038/srep12276',
        pmcid='PMC4505325',
    ))


def test_caching_works_as_expected():
    fetch.summary.cache_clear()
    assert fetch.summary.cache_info().hits == 0
    assert fetch.summary.cache_info().misses == 0
    assert lookup('PMID:375006')
    assert fetch.summary.cache_info().hits == 0
    assert fetch.summary.cache_info().misses == 1
    for count in range(10):
        assert lookup('PMID:375006')
        assert fetch.summary.cache_info().hits == count + 1
        assert fetch.summary.cache_info().misses == 1


def test_can_find_unidexable_publication():
    ref = lookup("1903816")
    assert ref == attr.asdict(Reference(
        authors='Leung J, Sinclair DA, Hayashi S, Tener GM, Grigliatti TA.',
        location='J Mol Biol 219(2):175-188 (1991)',
        title=(
            'Informational redundancy of tRNA(4Ser) and tRNA(7Ser) genes in '
            'Drosophila melanogaster and evidence for intergenic recombination'
        ),
        pmid=1903816,
        doi='10.1016/0022-2836(91)90560-s',
        pmcid=None,
    ))
