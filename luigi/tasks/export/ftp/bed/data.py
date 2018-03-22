# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from tasks.config import db
from tasks.config import export

from rnacentral.export.ftp import bed

class ChromSizes(luigi.Task):
    taxid = luigi.Parameter(default=9606)

    def output(self):
        return luigi.LocalTarget(export().bed('%s.chrom.sizes' % self.taxid))

    def run(self):
        bed.get_chrom_sizes(db(), taxid=self.taxid, output=self.output().path)


class BedToBigBed(luigi.Task):
    taxid = luigi.Parameter(default=9606)

    def requires(self):
        return [BedFile(taxid=self.taxid), ChromSizes(taxid=self.taxid)]

    def output(self):
        return luigi.LocalTarget(export().bed('rnacentral-%i.bigbed' % self.taxid))

    def run(self):
        bed.convert_to_bigbed(BedFile(taxid=self.taxid).output().path, ChromSizes(taxid=self.taxid).output().path, self.output().path)


class BedFile(luigi.Task):
    taxid = luigi.Parameter(default=9606)

    def requires(self):
        return BedDataDump(taxid=self.taxid)

    def output(self):
        return luigi.LocalTarget(export().bed('rnacentral-%i.bed' % self.taxid))

    def run(self):
        with self.input().open('r') as raw, self.output().open('w') as out:
            bed.make_bed_file(raw, out)
        bed.sort_bed_file(self.output().path)


class BedDataDump(luigi.Task):
    taxid = luigi.Parameter(default=9606)

    def output(self):
        return luigi.LocalTarget(export().bed('rnacentral-%i.tsv' % self.taxid))

    def run(self):
        with self.output().open('w') as raw:
            bed.coordinates(db(), raw, taxid=self.taxid)


class BedWrapper(luigi.Task):

    def requires(self):
        tasks = []
        for taxid in bed.get_taxids_with_genomic_mapping():
            tasks.append(BedToBigBed(taxid=taxid))
        return tasks
