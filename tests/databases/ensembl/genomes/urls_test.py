# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import typing as ty
from functools import lru_cache

import pytest

from rnacentral_pipeline.databases.ensembl.data import Division
from rnacentral_pipeline.databases.ensembl.genomes import urls


@lru_cache()
def species(division: Division) -> ty.List[str]:
    found = urls.urls_for(division, 'ftp.ensemblgenomes.org')
    return [f[1] for f in found]


@pytest.mark.parametrize('division,expected,found', [
    (Division.fungi, 'aspergillus_oryzae_gca_002007945', False),
    (Division.fungi, 'aspergillus_fumigatus_var_rp_2014_gca_000731615', False),
    (Division.fungi, 'aspergillus_oryzae', True),
])
def test_can_generate_reasonable_species(division, expected, found):
    if found:
        assert expected in species(division)
    else:
        assert expected not in species(division)
