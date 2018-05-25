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

import os
import tempfile

import luigi
import pytest
import requests

from tasks.utils.http import download



def test_can_download_file():
    filename = tempfile.mktemp()
    assert not os.path.exists(filename)
    download('http://rnacentral.org', filename)
    assert os.path.exists(filename)
    with open(filename, 'r') as out:
        assert next(out) == '<!--\n'


def test_can_download_file_to_local_target():
    filename = tempfile.mktemp()
    assert not os.path.exists(filename)
    output = luigi.LocalTarget(filename)
    download('http://rnacentral.org', output)
    assert os.path.exists(filename)
    with open(filename, 'r') as out:
        assert next(out) == '<!--\n'


def test_writes_nothing_if_download_fails():
    filename = tempfile.mktemp()
    assert not os.path.exists(filename)
    with pytest.raises(requests.HTTPError):
        download('http://rnacentral.org/nonsense-url', filename)
    assert not os.path.exists(filename)
