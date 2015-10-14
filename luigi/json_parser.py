"""
Copyright [2009-2014] EMBL-European Bioinformatics Institute
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

import json
import luigi
import os

from csv_output_writer import CsvOutputWriter


class JsonParser(luigi.Task, CsvOutputWriter):  # pylint: disable=W0232
    """
    Base class for all json processors.
    """
    input_file = luigi.Parameter()
    destination = luigi.Parameter(default='/tmp')
    test = luigi.BoolParameter(default=False, significant=False)
    data = [] # json input data
    entries = [] # RNAcentralEntry objects

    def get_accession(self, entry, database):
        """
        Generate unique accession.
        Example: CM000670.1:104203457..104215098:ncRNA:lncipedia
        """
        return ('{parent_accession}.{seq_version}:{start}..{stop}:{ncRNA_type}:'
                '{database}:{external_id}').format(
            parent_accession=entry.parent_accession,
            seq_version=entry.seq_version,
            start=entry.feature_location_start,
            stop=entry.feature_location_end,
            ncRNA_type=entry.feature_type,
            database=database,
            external_id=entry.primary_id
        )

    def get_description(self, entry):
        """
        Get entry description.
        """
        return '{species} {product}'.format(
            species=entry.species,
            product=entry.product
        )

    def get_feature_start_end(self, assembly_info):
        """
        Get feature start and end based on all exons.
        """
        primary_start = primary_end = []
        for exon in assembly_info:
            primary_start.append(exon['primary_start'])
            primary_end.append(exon['primary_end'])
        return (min(primary_start), max(primary_end))

    def parse_input_data(self):
        """
        Parse the input json file into a Python object.
        """
        with self.input().open('r') as f: # pylint: disable=E1101
            self.data = json.load(f)

    def write_output(self):
        """
        Create output files.
        """
        self.format_sequence_data(self.entries, self.output()['short'].fn, self.output()['long'].fn) # pylint: disable=E1101
        self.format_references(self.entries, self.output()['refs'].fn) # pylint: disable=E1101
        self.format_accession_info(self.entries, self.output()['ac_info'].fn) # pylint: disable=E1101
        self.format_genomic_locations(self.entries, self.output()['genomic_locations'].fn) # pylint: disable=E1101

    def requires(self):
        """
        Make sure that the input file exists.
        """
        filename = self.input_file
        class InputFile(luigi.Task):  # pylint: disable=W0232
            """
            Wrap input file in a luigi task.
            """
            def output(self):
                """
                """
                return luigi.LocalTarget(filename)
        return InputFile()

    def run(self):
        """
        Main task code.
        """
        self.parse_input_data()
        self.create_rnacentral_entries() # pylint: disable=E1101
        self.write_output()

    def output(self):
        """
        Generate output file paths and return them as a dictionary.
        """
        if not os.path.exists(self.destination):
            os.makedirs(self.destination)
        outputs = {}
        folders = ['short', 'long', 'ac_info', 'genomic_locations', 'refs']
        for folder in folders:
            subpath = os.path.join(self.destination, folder)
            if not os.path.exists(subpath):
                os.makedirs(subpath)
            filename = ''.join([self.database, '_', os.path.basename(self.input().fn).split('.')[0], '_', folder, '.csv']) # pylint: disable=E1101
            outputs[folder] = luigi.LocalTarget(os.path.join(subpath, filename))
        return outputs
