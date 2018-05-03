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
import shutil
import tempfile
from contextlib import contextmanager

import pytest

from tasks.utils.fetch import fetch


@contextmanager
def tmp_dir():
    dirname = tempfile.mkdtemp()
    try:
        yield dirname
    finally:
        shutil.rmtree(dirname)



@pytest.mark.parametrize('filename', [
    'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/_README.TXT',
    'http://www.google.com',
    'readme.md',
])
def test_can_copy_a_over(filename):
    with tempfile.NamedTemporaryFile() as tmp:
        fetch(filename, tmp.name)
        tmp.flush()
        with open(tmp.name, 'r') as final:
            assert final.read()


def test_can_copy_a_file_into_a_directory():
    with tmp_dir() as tdir:
        with tempfile.NamedTemporaryFile() as tmp:
            tmp.write('Hello\n')
            tmp.flush()
            fetch(tmp.name, tdir)
            final = os.path.join(tdir, os.path.basename(tmp.name))
            assert os.path.exists(final)
            with open(final, 'r') as raw:
                assert raw.read() == 'Hello\n'
            assert os.path.exists(tmp.name)


@pytest.mark.xfail(reason='Not sure what the correct test/behavior is')
def test_can_copy_a_directory():
    with tmp_dir() as tdir:
        os.system('touch %s' % (os.path.join(tdir, 't')))
        os.makedirs(os.path.join(tdir, 'a', 'b'))
        with tmp_dir() as target_dir:
            fetch(tdir, target_dir)
            assert os.path.exists(os.path.join(target_dir, 't'))
            assert os.path.exists(os.path.join(target_dir, 'a', 'b'))


@pytest.mark.skip
@pytest.mark.parametrize('location', [
    'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/not-here.TXT'
    'a-path-that-does-not-exist/yup',
])
def test_complains_given_bad_locations(location):
    filename = 'should-never-exist'
    with pytest.raises(Exception):
        fetch(location, filename)

    assert not os.path.exists(filename)
