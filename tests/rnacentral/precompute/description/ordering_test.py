# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.data import Database
from rnacentral_pipeline.rnacentral.precompute.description.species_specific import (
    ORDERING,
)


# A sequence only a missing database provides gets no description, failing the run
def test_every_database_can_provide_a_description():
    assert [d.name for d in Database if d not in ORDERING] == []
