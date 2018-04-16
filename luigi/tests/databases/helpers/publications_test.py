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

from databases.data import Reference
import databases.helpers.publications as pub


def test_can_fetch_publication():
    assert pub.summary(28815543) == {
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
        "pubType": "review; journal article; research support, n.i.h., extramural; ",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "citedByCount": 0,
        "hasReferences": "N",
        "hasTextMinedTerms": "N",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2017-01-01"
    }


def test_complains_given_bad_pmid():
    with pytest.raises(pub.UnknownPmid):
        pub.summary(-1)


def test_can_build_reference():
    assert attr.asdict(pub.reference('a', 27858507)) == attr.asdict(Reference(
        accession='a',
        authors=u"Lünse CE, Weinberg Z, Breaker RR.",
        location='RNA Biol 14(11):1499-1507 (2017)',
        title=(
            'Numerous small hammerhead ribozyme variants associated with '
            'Penelope-like retrotransposons cleave RNA as dimers'
        ),
        pmid=27858507,
        doi='10.1080/15476286.2016.1251002'
    ))


def test_can_deal_with_unicode():
    reference = pub.reference('a', 27334534)
    assert u'\xa0' not in reference.title
    assert reference.md5() == 'a84bed065b6f62d0c096d8bd7547b578'


def test_builds_correction_location():
    assert attr.asdict(pub.reference('bob', 26184978)) == attr.asdict(Reference(
        accession='bob',
        authors='Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.',
        location='Sci Rep 5:12276 (2015)',
        title=(
            'MiR-135b-5p and MiR-499a-3p Promote Cell '
            'Proliferation and Migration in Atherosclerosis by Directly '
            'Targeting MEF2C'
        ),
        pmid=26184978,
        doi='10.1038/srep12276',
    ))


def test_can_handle_missing_volume():
    assert attr.asdict(pub.reference('other', 27389411)) == attr.asdict(Reference(
        accession='other',
        authors='Carr G, Barrese V, Stott JB, Povstyan OV, Jepps TA, Figueiredo HB, Zheng D, Jamshidi Y, Greenwood IA.',
        location='Cardiovasc Res (2016)',
        title='MicroRNA-153 targeting of KCNQ4 contributes to vascular dysfunction in hypertension',
        pmid=27389411,
        doi='10.1093/cvr/cvw177',
    ))
