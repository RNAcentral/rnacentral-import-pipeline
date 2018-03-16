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

import luigi

from tasks.config import export
from tasks.utils.files import atomic_output

README = """
===================================================================
RNAcentral MD5 Data
===================================================================

Each RNAcentral id corresponds to an md5 hash calculated
based on the uppercase DNA version of that sequence.

* md5.tsv.gz
Tab-separated file with RNAcentral ids and md5 hashes
of their corresponding sequences.
This file can be used to look up RNAcentral ids for a set of sequences
given their md5.

* example.txt
A small file showing the first few md5 entries.
"""


class Readme(luigi.Task):
    def output(self):
        return luigi.LocalTarget(export().md5('readme.txt'))

    def run(self):
        with atomic_output(self.output()) as out:
            out.write(README)
