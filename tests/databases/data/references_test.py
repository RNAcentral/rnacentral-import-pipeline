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
import requests

from rnacentral_pipeline.databases import data


def test_reference_can_handle_non_unicode():
    ref = data.Reference(
        authors='Ahn S., Jin T.E., Chang D.H., Rhee M.S., Kim H.J., Lee S.J., Park D.S., Kim B.C.',
        location='Int. J. Syst. Evol. Microbiol. 66(9):3656-3661(2016).',
        title='Agathobaculum butyriciproducens gen. nov. \xa0sp. nov., a strict anaerobic, butyrate-producing gut bacterium isolated from human faeces and reclassification of Eubacterium desmolans as Agathobaculum desmolans comb. nov',
        pmid=27334534,
        doi=None,
        pmcid=None,
    )
    assert ref.md5() == '1c1aa1c716a1ae7fd6ba0747d3e166e0'


@pytest.mark.parametrize('pmid', [
    'PMID:30715521',
    '30715521',
    '  30715521',
    30715521,
])
def test_can_build_id_reference_for_simple_pmids(pmid):
    assert attr.asdict(data.IdReference.build(pmid)) == attr.asdict(data.IdReference(
        namespace=data.KnownServices.pmid,
        external_id='30715521',
    ))


@pytest.mark.parametrize('pmc', [
    'PMC4648659',
    'pmc4648659',
    'PMCID:PMC4648659',
    'pmcid:pmc4648659',
])
def test_can_build_id_reference_for_pmcid(pmc):
    assert attr.asdict(data.IdReference.build(pmc)) == attr.asdict(data.IdReference(
        namespace=data.KnownServices.pmcid,
        external_id='PMC4648659'
    ))


@pytest.mark.parametrize('raw_data,title', [
    ('PMC1172314', 'Books also Received.'),
    ('PMC3131748', 'Abstracts - Invited Speakers.'),
    ('PMC3131791', 'Abstracts - Poster Presentations.'),
    ('PMID:30715521', 'LncBook: a curated knowledgebase of human long non-coding RNAs.'),
    ('doi:10.1007/s11524-007-9234-y', 'Section I: Oral Sessions.'),
    ('doi:10.1016/S1470-2045(16)30240-6', "Genome-wide association studies in oesophageal adenocarcinoma and Barrett's oesophagus: a large-scale meta-analysis."),
    ('doi:10.1016/S1672-0229(04)02021-2', "Application of proteomics in the study of tumor metastasis."),
    ('doi:10.1083/jcb.1851iti3', 'What do kidneys and embryonic fish skin have in common?'),
    ('doi:10.1083/jcb.2111iti1', 'mTORC2 tips the balance in cell survival.'),
    ('doi:10.1093/database/baw138', 'IRNdb: the database of immunologically relevant non-coding RNAs.'),
    ('doi:10.1093/ofid/ofx163.867', 'The Expression of hsp-miRNA-200b-3p and -200c-3p in Human Cytomegalovirus-infected Formalin-Fixed, Paraffin-Embedded Tissues.'),
    ('doi:10.1093/schbul/sby014.092', '23.2 NETRIN-1 RECEPTORS CONTROL MESOCORTICAL DOPAMINE CONNECTIVITY IN ADOLESCENCE.'),
    ('doi:10.1097/01.WOX.0000411770.14047.89', '25 Role of Myeloid Derived Suppressor Cells in Asthma.'),
    ('doi:10.1097/MD.0000000000002371', 'UMIB Summit 2015.'),
    ('doi:10.1111/cas.12358', 'In This Issue.'),
    ('doi:10.1111/cas.12385', 'In This Issue.'),
    ('doi:10.1177/2050640615601611', 'UEG Week 2015 Oral Presentations.'),
    ('doi:10.1371/journal.pbio.0020114', 'Exploring Small RNA Function.'),
    ('doi:10.1534/g3.116.036848', 'Meeting Report: The Allied Genetics Conference 2016.'),
    ('pmcid:PMC1480519', 'Fellowships, Grants, & Awards.'),
    ('pmcid:PMC5064671', 'Abstracts - USICON 2016.'),
    ('pmid:17254355', 'Retroviral activation of the mir-106a microRNA cistron in T lymphoma.'),
    ('pmid:24531370', 'New insights into the promoterless transcription of DNA coligo templates by RNA polymerase III.'),
])
def test_can_query_for_expected_data(raw_data, title):
    ref = data.IdReference.build(raw_data)
    response = requests.get(ref.external_url())
    assert response.json()['hitCount'] == 1
    assert response.json()['resultList']['result'][0]['title'] == title
