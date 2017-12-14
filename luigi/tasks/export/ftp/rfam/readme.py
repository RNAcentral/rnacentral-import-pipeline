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
from luigi.local_target import atomic_file

from tasks.config import export

README = """
===================================================================
RNAcentral Rfam Annotation Data
===================================================================

This directory contains mappings between RNAcentral sequences and Rfam
annotations for those sequences.

* rfam_annotations.tsv.gz
    A tab separated file of all active RNAcentral sequences and their
    associated Rfam annotations. The file has the format:

    URS-Id Rfam-Model-Id Score E-value Sequence-Start Sequence-Stop Model-Start Model-Stop Rfam-Model-Description

    The start and stop positions are 0 indexed.

* example.txt
    A small example file showing the format of rfam_annotations.tsv.gz.
"""


class RfamAnnotationsReadMe(luigi.Task):
    def output(self):
        return luigi.LocalTarget(export().ftp('rfam', 'readme.txt'))

    def run(self):
        with atomic_file(self.output().fn) as out:
            out.write(README)
