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

import tempfile

import six

from rnacentral_pipeline import utils


def test_can_serialize_stream():
    data = ["a", 1, 2, 3, 4]
    with tempfile.NamedTemporaryFile() as tmp:
        utils.pickle_stream(data, tmp)
        tmp.seek(0)
        result = utils.unpickle_stream(tmp)
        assert not isinstance(result, list)
        for index, obj in enumerate(result):
            assert data[index] == obj
