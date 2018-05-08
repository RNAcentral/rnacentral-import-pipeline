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

import csv
import json
import os

from copy import copy

import luigi

from tasks.config import db
from tasks.config import export
from tasks.genome_mapping.genome_mapping_tasks import get_mapped_assemblies
from tasks.genome_mapping.genome_mapping_tasks import GetAssemblyInfoJson
from rnacentral.export.ftp import bed


class ChromSizesUCSC(luigi.Task):
    """
    Retrieve a chromosome sizes file for a genome from UCSC.
    """
    assembly_ucsc = luigi.Parameter()

    def output(self):
        filename = '{}.chrom.sizes'.format(self.assembly_ucsc)
        return luigi.LocalTarget(export().bed(filename))

    def run(self):
        bed.get_chrom_sizes(db(), assembly_ucsc=self.assembly_ucsc, output=self.output().path)


class ChromSizesEnsembl(luigi.Task):
    """
    Generate a chromosome sizes file locally based on Ensembl data.
    """
    species = luigi.Parameter()
    division = luigi.Parameter()

    def output(self):
        filename = '{}.chrom.sizes'.format(self.species.capitalize())
        return luigi.LocalTarget(export().bed(filename))

    def requires(self):
        return GetAssemblyInfoJson(species=self.species, division=self.division)

    def run(self):
        data = json.load(self.input().open('r'))
        with self.output().open('w') as outfile:
            for region in data['top_level_region']:
                name = bed.format_chromosome_name(region['name'], region)
                line = '{}\t{}\n'.format(name, region['length'])
                outfile.write(line)


class BedDataDump(luigi.Task):
    """
    Export a tsv file with genome coordinates that will be processed into bed
    format.
    """
    assembly_id = luigi.Parameter()
    species = luigi.Parameter()
    taxid = luigi.IntParameter()

    def output(self):
        filename = '{}.{}.tsv'.format(self.species.capitalize(), self.assembly_id)
        return luigi.LocalTarget(export().bed(filename))

    def run(self):
        with self.output().open('w') as raw:
            bed.export_blat_coordinates(db(), raw, taxid=self.taxid)
            bed.export_ensembl_coordinates(db(), raw, taxid=self.taxid)


class BedFile(luigi.Task):
    """
    Generate a bed file.
    """
    assembly_id = luigi.Parameter()
    species = luigi.Parameter()
    taxid = luigi.IntParameter()
    division = luigi.Parameter()

    def requires(self):
        return {
            'data': BedDataDump(taxid=self.taxid,
                                assembly_id=self.assembly_id,
                                species=self.species),
            'chrom_sizes': ChromSizesEnsembl(species=self.species,
                                             division=self.division),
            'ensembl_info': GetAssemblyInfoJson(species=self.species,
                                                division=self.division)
        }

    def output(self):
        filename = '{}.{}.bed'.format(self.species.capitalize(), self.assembly_id)
        return luigi.LocalTarget(export().bed(filename))

    def get_ensembl_json_data(self):
        """
        Restructure Ensembl Json object to speed up look up by making
        `top_level_region` a dictionary instead of an array.
        """
        data = json.load(self.input()['ensembl_info'].open('r'))
        new_data = copy(data)
        new_data['top_level_region'] = {}
        for region in data['top_level_region']:
            new_data['top_level_region'][region['name']] = region
        return new_data

    def run(self):
        with self.input()['data'].open('r') as raw, self.output().open('w') as out:
            bed.make_bed_file(raw, out, self.get_ensembl_json_data())
        bed.sort_bed_file(self.output().path)


class BedToBigBed(luigi.Task):
    """
    Convert bed to bigbed format.
    """
    assembly_ucsc = luigi.Parameter()
    assembly_id = luigi.Parameter()
    species = luigi.Parameter()
    taxid = luigi.IntParameter()
    division = luigi.Parameter()

    def requires(self):
        return {
            'bed': BedFile(taxid=self.taxid,
                           assembly_id=self.assembly_id,
                           division=self.division,
                           species=self.species),
            'chromsizes': ChromSizesEnsembl(species=self.species,
                                            division=self.division),
        }

    def output(self):
        filename = '{}.{}.bigBed'.format(self.species.capitalize(), self.assembly_id)
        return luigi.LocalTarget(export().bed(filename))

    def run(self):
        bed.convert_to_bigbed(self.input()['bed'].path,
                              self.input()['chromsizes'].path,
                              self.output().path)


class BigBedWrapper(luigi.WrapperTask):
    """
    Generate all bigbed files.
    """
    species = luigi.Parameter(default='')

    def requires(self):
        for assembly in get_mapped_assemblies():
            if self.species and self.species != assembly['species']:
                continue
            yield BedToBigBed(taxid=assembly['taxid'],
                              assembly_id=assembly['assembly_id'],
                              assembly_ucsc=assembly['assembly_ucsc'],
                              division=assembly['division'],
                              species=assembly['species'])


class BedWrapper(luigi.WrapperTask):
    """
    Generate all bed files.
    """
    species = luigi.Parameter(default='')

    def requires(self):
        for assembly in get_mapped_assemblies():
            if self.species and self.species != assembly['species']:
                continue
            yield BedFile(taxid=assembly['taxid'],
                          assembly_id=assembly['assembly_id'],
                          division=assembly['division'],
                          species=assembly['species'])


class BedAndBigBedWrapper(luigi.WrapperTask):
    """
    A task for generating all bed and bigbed files.
    """
    def requires(self):
        yield BigBedWrapper()
        yield BedWrapper()
