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
import string
import random
import tempfile

import luigi

from tasks.utils.files import atomic_output


def test_writes_file():
    filename = tempfile.mktemp()
    assert not os.path.exists(filename)
    with atomic_output(filename) as out:
        out.write("hi")

    assert os.path.exists(filename)
    with open(filename, 'rb') as data:
        assert data.read() == "hi"
    os.remove(filename)


def test_creates_directory():
    filename = tempfile.mktemp()
    dirs = ''.join(random.choice(string.ascii_lowercase) for i in xrange(10))
    filename = os.path.join(os.path.dirname(filename), dirs, 't')
    dirname = os.path.dirname(filename)
    assert not os.path.exists(dirname)
    with atomic_output(filename) as out:
        out.write("hi")
    assert os.path.exists(dirname)
    with open(filename, 'rb') as data:
        assert data.read() == "hi"
    os.remove(filename)
    os.rmdir(dirname)


def test_will_not_write_if_exception():
    filename = tempfile.mktemp()
    assert not os.path.exists(filename)
    try:
        with atomic_output(filename) as out:
            out.write("hi")
            raise Exception("Failing")
    except:
        pass

    assert not os.path.exists(filename)


def test_can_work_with_local_target():
    filename = tempfile.mktemp()
    assert not os.path.exists(filename)
    output = luigi.LocalTarget(filename)
    with atomic_output(output) as out:
        out.write("hi")

    assert os.path.exists(filename)
    with open(filename, 'rb') as data:
        assert data.read() == "hi"
    os.remove(filename)
