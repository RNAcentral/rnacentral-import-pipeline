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

import luigi
from luigi import LocalTarget
from luigi.local_target import atomic_file

import requests

from tasks.config import output

URL = 'http://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt'


class MgiDownload(luigi.Task):  # pylint: disable=R0904
    """
    This will download the MGI MRK_Sequence.rpt file for future processing.
    """

    def output(self):
        path = os.path.join(output().base, 'mgi', 'MRK_Sequence.rpt')
        return LocalTarget(path)

    def run(self):
        filename = self.output().fn
        dirname = os.path.dirname(filename)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        response = requests.get(URL)
        response.raise_for_status()
        with atomic_file(filename) as out:
            out.write(response.text)
