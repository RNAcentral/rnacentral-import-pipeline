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
from datetime import date
import xml.etree.ElementTree as ET

import luigi

from rnacentral.utils import upi_ranges

from tasks.config import db
from tasks.config import export
from tasks.config import output

from .validate import ValidateAndCompressSearchChunk

NOTE = """
release=8.0
release_date={date}
entries={count}
"""


class SearchReleaseNote(luigi.Task):
    def requires(self):
        config = export()
        for start, stop in upi_ranges(db(), config.search_export_size):
            yield ValidateAndCompressSearchChunk(min=start, max=stop)

    def output(self):
        config = output()
        filename = os.path.join(config.search_files, 'release_note.txt')
        return luigi.LocalTarget(filename)

    def counts(self):
        total = 0
        for requirement in self.requires():
            with requirement.output().open('r') as raw:
                tree = ET.parse(raw)
                root = tree.getroot()
                count = root.find('entry_count')
                if not count or not count.text:
                    raise ValueError("Could not find a count")
                if int(count.text) <= 0:
                    raise ValueError("Count must be positive")
                total += int(count.text)
        return total

    def note(self):
        return NOTE.format(
            date=date.strftime('%d-%B-%Y'),
            count=self.counts(),
        )

    def run(self):
        with self.output().open('w') as raw:
            raw.write(self.note())
