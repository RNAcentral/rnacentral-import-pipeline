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
import shutil
from contextlib import contextmanager

from unittest import TestCase

from tasks import config as conf


def touch(*parts):
    path = os.path.join(*parts)
    with open(path, 'a'):
        os.utime(path, None)


class EnaTest(TestCase):
    @contextmanager
    def tmp_base_config(self):
        config = conf.ena()
        tmpdir = tempfile.mkdtemp()
        config.base = tmpdir
        try:
            yield config
        finally:
            shutil.rmtree(tmpdir)

    def test_raw_ncr_files_does_not_include_empty_directories(self):
        with self.tmp_base_config() as config:
            os.makedirs(os.path.join(config.base, 'wgs', 'aa'))
            os.makedirs(os.path.join(config.base, 'wgs', 'ab'))
            touch(config.base, 'wgs', 'aa', 'a.ncr.gz')

            assert config.raw_ncr_files() == [
                os.path.join(config.base, 'wgs', 'aa')
            ]

    def test_raw_ncr_files_does_not_include_fasta_directories(self):
        with self.tmp_base_config() as config:
            os.makedirs(os.path.join(config.base, 'std', 'fasta'))
            os.makedirs(os.path.join(config.base, 'rel', 'fasta'))
            touch(config.base, 'std', 'fasta', 'ex.fasta.gz')
            touch(config.base, 'std', 'std_aa.ncr.gz')
            touch(config.base, 'rel', 'rel_1.ncr.gz')

            assert sorted(config.raw_ncr_files()) == sorted([
                os.path.join(config.base, 'std', 'std_aa.ncr.gz'),
                os.path.join(config.base, 'rel', 'rel_1.ncr.gz'),
            ])
