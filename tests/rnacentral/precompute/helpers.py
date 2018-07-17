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


import os

from rnacentral_pipeline.rnacentral.precompute.process import as_sequences
from rnacentral_pipeline.rnacentral.precompute.data import Sequence

from tests.helpers import run_range_as_single


def load_data(rna_id):
    path = os.path.join('files', 'precompute', 'query.sql')
    data = run_range_as_single(rna_id, path)
    return next(as_sequences([data]))
