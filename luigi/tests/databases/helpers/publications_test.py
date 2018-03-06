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
    assert '\xa0' not in reference.title
    assert reference.md5() == 'a84bed065b6f62d0c096d8bd7547b578'
