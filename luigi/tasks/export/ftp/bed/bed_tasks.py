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
from tasks.genome_mapping.genome_mapping_tasks import get_mapped_assemblies
from rnacentral.export.ftp import bed


class ChromSizes(luigi.Task):
    assembly_ucsc = luigi.Parameter(default='hg38')

    def output(self):
        filename = '{}.chrom.sizes'.format(self.assembly_ucsc)
        return luigi.LocalTarget(export().bed(filename))

    def run(self):
        bed.get_chrom_sizes(db(), assembly_ucsc=self.assembly_ucsc, output=self.output().path)


class BedDataDump(luigi.Task):
    assembly_id = luigi.Parameter(default='GRCh38')
    species = luigi.Parameter(default='homo_sapiens')
    taxid = luigi.IntParameter(default=9606)

    def output(self):
        filename = '{}.{}.tsv'.format(self.species.capitalize(), self.assembly_id)
        return luigi.LocalTarget(export().bed(filename))

    def run(self):
        with self.output().open('w') as raw:
            bed.export_blat_coordinates(db(), raw, taxid=self.taxid)
            bed.export_ensembl_coordinates(db(), raw, taxid=self.taxid)


class BedFile(luigi.Task):
    assembly_id = luigi.Parameter(default='GRCh38')
    species = luigi.Parameter(default='homo_sapiens')
    taxid = luigi.IntParameter(default=9606)

    def requires(self):
        return BedDataDump(taxid=self.taxid,
                           assembly_id=self.assembly_id,
                           species=self.species)

    def output(self):
        filename = '{}.{}.bed'.format(self.species.capitalize(), self.assembly_id)
        return luigi.LocalTarget(export().bed(filename))

    def run(self):
        with self.input().open('r') as raw, self.output().open('w') as out:
            bed.make_bed_file(raw, out)
        bed.sort_bed_file(self.output().path)


class BedToBigBed(luigi.Task):
    assembly_ucsc = luigi.Parameter(default='hg38')
    assembly_id = luigi.Parameter(default='GRCh38')
    species = luigi.Parameter(default='homo_sapiens')
    taxid = luigi.IntParameter(default=9606)

    def requires(self):
        return {
            'bed': BedFile(taxid=self.taxid,
                           assembly_id=self.assembly_id,
                           species=self.species),
            'chromsizes': ChromSizes(assembly_ucsc=self.assembly_ucsc),
        }

    def output(self):
        filename = '{}.bigBed'.format(self.assembly_ucsc)
        return luigi.LocalTarget(export().bed(filename))

    def run(self):
        bed.convert_to_bigbed(self.input()['bed'].path,
                              self.input()['chromsizes'].path,
                              self.output().path)


class BigBedWrapper(luigi.WrapperTask):
    def requires(self):
        for assembly in get_mapped_assemblies():
            if assembly['assembly_ucsc']:
                yield BedToBigBed(taxid=assembly['taxid'],
                                  assembly_id=assembly['assembly_id'],
                                  assembly_ucsc=assembly['assembly_ucsc'],
                                  species=assembly['species'])


class BedWrapper(luigi.WrapperTask):
    def requires(self):
        for assembly in get_mapped_assemblies():
            yield BedFile(taxid=assembly['taxid'],
                          assembly_id=assembly['assembly_id'],
                          species=assembly['species'])


class BedAndBigBedWrapper(luigi.WrapperTask):
    def requires(self):
        yield BigBedWrapper()
        yield BedWrapper()
