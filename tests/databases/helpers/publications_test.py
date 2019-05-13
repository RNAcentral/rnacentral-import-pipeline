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

import csv
import sqlite3
import tempfile

import six
import attr
import pytest
import shutil

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.data import IdReference
import rnacentral_pipeline.databases.helpers.publications as pub


def lookup(ref_id):
    return attr.asdict(pub.query_pmc(IdReference.build(ref_id)))


@pytest.fixture(scope='module')
def indexed_db_path():
    tmp = tempfile.mkdtemp()
    pub.index_xml_directory('data/publications/', tmp)
    yield tmp
    shutil.rmtree(tmp)


@pytest.fixture(scope='module')
def indexed(indexed_db_path):
    pub.index_xml_directory('data/publications/', indexed_db_path)
    with pub.Cache.build(indexed_db_path).open('r') as db:
        yield db


@pytest.mark.parametrize('raw_id', [
    28815543,
    '28815543',
    'PMC5890441',
    'doi:10.1007/978-981-10-5203-3_9',
])
def test_can_fetch_publication(raw_id):
    assert pub.summary(IdReference.build(raw_id)) == {
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
        "citedByCount": 3,
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
        pub.summary(IdReference.build(-1))


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
    reference = lookup(27334534)
    assert u'\xa0' not in reference['title']
    # assert reference.md5() == 'a84bed065b6f62d0c096d8bd7547b578'


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


@pytest.mark.parametrize('raw_id,namespace,external_id', [
    (26184978, 'pmid', '26184978'),
    ('pmid:26184978', 'pmid', '26184978'),
    ('PMID:26184978', 'pmid', '26184978'),
    ('doi:10.1038/srep12276', 'doi', '10.1038/srep12276'),
    ('DOI:10.1038/srep12276', 'doi', '10.1038/srep12276'),
    ('PMCID:PMC5785218', 'pmcid', 'PMC5785218'),
    ('PMC5785218', 'pmcid', 'PMC5785218'),
    ('pmcid:pmc5785218', 'pmcid', 'PMC5785218'),
    ('pmc5785218', 'pmcid', 'PMC5785218'),
])
def test_reference_builds_a_id_ref(raw_id, namespace, external_id):
    assert pub.reference(raw_id) == IdReference(
        namespace=namespace,
        external_id=external_id,
    )


def test_can_parse_xml_data_correctly():
    with open('data/publications/example.xml', 'r') as raw:
        data = list(pub.parse_xml(raw))
        assert len(data) == 1
        assert attr.asdict(data[0]) == attr.asdict(Reference(
            authors='Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.',
            location='Scientific reports 5:12276 (2015)',
            title=(
                'MiR-135b-5p and MiR-499a-3p Promote Cell '
                'Proliferation and Migration in Atherosclerosis by Directly '
                'Targeting MEF2C'
            ),
            pmid=26184978,
            doi='10.1038/srep12276',
            pmcid='PMC4505325',
        ))


@pytest.mark.parametrize('raw_id', [
    ('26184978'),
    ('DOI:10.1038/srep12276'),
    ('PMC4505325'),
    ('PMCID:PMC4505325'),
    (26184978),
])
def test_can_query_indexed_data_correctly(indexed, raw_id):
    id_ref = pub.reference(raw_id)
    db = indexed[id_ref.namespace]
    assert attr.asdict(db.get(id_ref)) == attr.asdict(Reference(
            authors='Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.',
            location='Scientific reports 5:12276 (2015)',
            title=(
                'MiR-135b-5p and MiR-499a-3p Promote Cell '
                'Proliferation and Migration in Atherosclerosis by Directly '
                'Targeting MEF2C'
            ),
            pmid=26184978,
            doi='10.1038/srep12276',
            pmcid='PMC4505325',
        ))


@pytest.mark.parametrize('raw_id', [
    'doi:10.1007/bf00271669',
    'PMID:375006',
])
def test_can_query_with_fallback(indexed, raw_id):
    id_ref = pub.reference(raw_id)
    db = indexed[id_ref.namespace]
    ref = db.get(id_ref, allow_fallback=True)
    assert attr.asdict(ref) == attr.asdict(Reference(
        authors='Macino G, Tzagoloff A.',
        location='Mol Gen Genet 169(2):183-188 (1979)',
        title='Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae',
        pmid=375006,
        doi='10.1007/bf00271669',
        pmcid=None,
    ))


@pytest.mark.parametrize('raw_id', [
    26184978,
    'pmid:26184978',
    'PMID:26184978',
    'doi:10.1038/srep12276',
    'DOI:10.1038/srep12276',
])
def test_can_write_using_specified_columns(indexed_db_path, raw_id):
    out = six.moves.StringIO()
    content = 'something,%s,other\n' % raw_id
    raw = six.moves.StringIO(content)
    pub.write_file_lookup(indexed_db_path, raw, out, column=1, allow_fallback=False)
    out.seek(0)
    assert list(csv.reader(out)) == [
        [
            "03ea261d3cd947fde1cc8328a4c08127", 
            "something", 
            "other", 
            'Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.', 
            'Scientific reports 5:12276 (2015)', 
            (
                'MiR-135b-5p and MiR-499a-3p Promote Cell '
                'Proliferation and Migration in Atherosclerosis by Directly '
                'Targeting MEF2C'
            ),
            '26184978',
            '10.1038/srep12276'
        ]
    ]


@pytest.mark.parametrize('raw_id', [
    '375006',
    'PMID:375006',
    'DOI:10.1007/bf00271669',
])
def test_can_write_using_specified_columns_and_allow_fallback(indexed_db_path, raw_id):

    with pytest.raises(pub.UnknownReference):
        with pub.Cache.build(indexed_db_path).open('r') as cache:
            id_ref = pub.reference(raw_id)
            db = cache[id_ref.namespace]
            ref = db.get(id_ref, allow_fallback=False)

    out = six.moves.StringIO()
    raw = six.moves.StringIO('something,%s,other\n' % raw_id)
    pub.write_file_lookup(indexed_db_path, raw, out, column=1, allow_fallback=True)
    out.seek(0)
    assert list(csv.reader(out)) == [
        [
            "da0c9805cab7efd21a0eed0e17a8223a", 
            "something", 
            "other", 
            'Macino G, Tzagoloff A.', 
            'Mol Gen Genet 169(2):183-188 (1979)', 
            'Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae',
            '375006',
            '10.1007/bf00271669'
        ]
    ]


def test_caching_works_as_expected():
    pub.summary.cache_clear()
    assert pub.summary.cache_info().hits == 0
    assert pub.summary.cache_info().misses == 0
    assert lookup('PMID:375006')
    assert pub.summary.cache_info().hits == 0
    assert pub.summary.cache_info().misses == 1
    for count in range(10):
        assert lookup('PMID:375006')
        assert pub.summary.cache_info().hits == count + 1
        assert pub.summary.cache_info().misses == 1
