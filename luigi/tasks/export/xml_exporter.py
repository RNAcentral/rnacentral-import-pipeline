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

from tasks.config import db
from tasks.config import output

from rnacentral.xml.exporter import export_range
from rnacentral.xml.exporter import as_document
from rnacentral.db import cursor


class XmlExporterTask(luigi.Task):  # pylint: disable=R0904
    """
    This is a task that will create an xml export for the given range of UPI's.
    """

    min = luigi.IntParameter()
    max = luigi.IntParameter()

    def output(self):
        config = output()
        filepattern = 'xml_export_{min}_{max}.xml'.format(
            min=self.min,
            max=self.max,
        )
        filename = os.path.join(config.search_files, filepattern)
        return LocalTarget(filename)

    def run(self):
        with atomic_file(self.output().fn) as raw:
            with cursor(db()) as cur:
                results = export_range(cur, self.min, self.max)
                as_document(results, raw)
