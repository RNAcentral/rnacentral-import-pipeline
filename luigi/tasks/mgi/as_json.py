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
import json

import luigi
from luigi import LocalTarget
from luigi.local_target import atomic_file

import attr

from tasks.mgi.download import MgiDownload
from databases.mgi.parser import rna_entries


class MgiToJson(luigi.Task):  # pylint: disable=R0904
    """
    This will convert the downloaded MGI data to a JSON file containing more
    structured data for the RNA entries.
    """

    def requires(self):
        return MgiDownload()

    def output(self):
        filename = self.requires().output().fn
        base = os.path.dirname(filename)
        return LocalTarget(os.path.join(base, 'rna.json'))

    def run(self):
        input_file = self.requires().output().fn
        data = []
        for entry in rna_entries(input_file):
            result = attr.asdict(entry)
            result['feature_type'] = entry.feature_type
            result['ncrna_class'] = entry.ncrna_class
            result['feature_location_start'] = entry.feature_location_start
            result['feature_location_end'] = entry.feature_location_end
            result['external_id'] = entry.accession
            for index, reference in enumerate(entry.references):
                result['references'][index]['md5'] = reference.md5()

            data.append(result)

        with atomic_file(self.output().fn) as out:
            json.dump(data, out)
