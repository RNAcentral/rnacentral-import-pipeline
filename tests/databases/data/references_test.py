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

from rnacentral_pipeline.databases import data


def test_reference_can_handle_non_unicode():
    ref = data.Reference(
        authors='Ahn S., Jin T.E., Chang D.H., Rhee M.S., Kim H.J., Lee S.J., Park D.S., Kim B.C.',
        location='Int. J. Syst. Evol. Microbiol. 66(9):3656-3661(2016).',
        title='Agathobaculum butyriciproducens gen. nov. \xa0sp. nov., a strict anaerobic, butyrate-producing gut bacterium isolated from human faeces and reclassification of Eubacterium desmolans as Agathobaculum desmolans comb. nov',
        pmid=27334534,
        doi=None,
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
        namespace='pmid',
        external_id='30715521',
    ))


def test_can_build_id_reference_for_pmcid():
    assert attr.asdict(data.IdReference.build('PMC4648659')) == attr.asdict(data.IdReference(
        namespace='pmc',
        external_id='PMC4648659'
    ))

