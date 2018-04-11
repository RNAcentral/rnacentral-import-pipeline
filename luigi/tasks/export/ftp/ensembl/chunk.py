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
from luigi import LocalTarget

from tasks.config import db
from tasks.config import export
from tasks.utils.files import atomic_output

from rnacentral.export.ftp import ensembl


class EnsemblExportChunk(luigi.Task):
    min = luigi.IntParameter()
    max = luigi.IntParameter()

    def output(self):
        config = export()
        return LocalTarget(
            config.ensembl_export(ensembl.range_filename(self.min, self.max))
        )

    def run(self):
        with atomic_output(self.output()) as raw:
            results = ensembl.range(db(), self.min, self.max)
            ensembl.write_and_validate(raw, results)
