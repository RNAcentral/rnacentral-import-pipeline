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
import gzip
import tarfile
import tempfile
import shutil
from contextlib import contextmanager

import pytest

from tasks.utils.compressed import expand


@pytest.fixture
def tarball():
    with tempfile.NamedTemporaryFile(suffix='.tar.gz', delete=False) as tmp:
        with tarfile.open(tmp.name, 'w:gz') as tar:
            tar.add('Dockerfile')
            tar.add('LICENSE')
        yield tmp.name


@contextmanager
def tmp_dir():
    tmp = tempfile.mkdtemp()
    yield tmp
    shutil.rmtree(tmp)


@pytest.fixture
def gzip_file():
    with tmp_dir() as dirname:
        filename = os.path.join(dirname, 'simple.json.gz')
        with gzip.open(filename, 'wb') as out:
            out.write('{"data": "useful"}')
        yield filename


def test_can_unpack_a_tar_file(tarball):
    with tmp_dir() as name:
        assert os.path.exists(os.path.join(name, 'Dockerfile')) is False
        assert os.path.exists(os.path.join(name, 'LICENSE')) is False
        print("TARGET " + name)
        expand(tarball, name)
        print(os.listdir(name))
        assert os.path.exists(os.path.join(name, 'Dockerfile'))
        assert os.path.exists(os.path.join(name, 'LICENSE'))


def test_returns_list_of_unpacked_tar_files(tarball):
    with tmp_dir() as name:
        assert list(expand(tarball, name)) == [
            os.path.join(name, 'Dockerfile'),
            os.path.join(name, 'LICENSE'),
        ]


def test_can_unpack_a_gzip_file(gzip_file):
    with tmp_dir() as name:
        target = os.path.join(name, 'simple.json')
        assert os.path.exists(target) is False
        assert expand(gzip_file, name)
        assert os.path.exists(target)


def test_returns_an_iterable_of_unzipped_files(gzip_file):
    with tmp_dir() as name:
        assert list(expand(gzip_file, name)) == [
            os.path.join(name, 'simple.json'),
        ]
