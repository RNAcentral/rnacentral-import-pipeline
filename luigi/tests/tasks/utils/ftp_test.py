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
import ftplib
import tempfile

import pytest

from tasks.utils import ftp


def test_can_download_file_from_ftp():
    path = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/_README.TXT'
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        ftp.download(path, tmp.name)

    with open(tmp.name) as output:
        contents = output.read()

    assert contents
    assert contents.startswith('\n#######')

    os.unlink(tmp.name)


def test_it_fails_given_unknown_path():
    path = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/not-here.TXT'
    with pytest.raises(ftplib.error_perm):
        with tempfile.NamedTemporaryFile() as tmp:
            ftp.download(path, tmp.name)


def test_it_does_not_create_file_if_failing():
    path = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/not-here.TXT'
    filename = 'should-not-exist'
    with pytest.raises(ftplib.error_perm):
        ftp.download(path, filename)

    assert not os.path.exists(filename)
