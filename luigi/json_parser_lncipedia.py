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

Usage:

python path/to/this/file.py JsonParserLncipedia
    --local-scheduler
    --destination /path/to/output/files
    --input-file /path/to/input/file.json

Optional parameters:
    --test # to process just the first few entries (default: False)
    --no_assembly_mapping # do not map LNCipedia data onto the latest assembly (default: False)
"""

import json
import luigi
import os
import requests
import time

from rnacentral_entry import RNAcentralEntry
from csv_output_writer import CsvOutputWriter


class JsonParserLncipedia(luigi.Task, CsvOutputWriter):  # pylint: disable=W0232
    """
    Luigi Task for converting LNCipedia Json file into csv files
    that can be loaded into the RNAcentral database.
    """
    input_file = luigi.Parameter()
    destination = luigi.Parameter(default='/tmp')
    test = luigi.BoolParameter(default=False, significant=False)
    no_assembly_mapping = luigi.BoolParameter(default=True, significant=False)
    data = [] # json input data
    entries = [] # RNAcentralEntry objects

    def get_accession(self, entry):
        """
        Generate unique accession.
        Example: CM000670.1:104203457..104215098:ncRNA:lncipedia
        """
        return '{parent_accession}.{seq_version}:{start}..{stop}:{ncRNA_type}:lncipedia'.format(
            parent_accession=entry.parent_accession,
            seq_version=entry.seq_version,
            start=entry.feature_location_start,
            stop=entry.feature_location_end,
            ncRNA_type=entry.feature_type
        )

    def get_description(self, entry):
        """
        Get entry description.
        """
        return '{species} {product}'.format(
            species=entry.species,
            product=entry.product
        )

    def map_genomic_assembly(self, rnacentral_entry):
        """
        Use Ensembl REST API to map LNCipedia entries from GRCh37 to GRCh38.
        """
        for exon in rnacentral_entry.assembly_info:
            lncipedia_assembly_id = 'GRCh37'
            ensembl_assembly_id = 'GRCh38'
            url = ('http://rest.ensembl.org/map/human/'
                   '{lncipedia_assembly_id}/{chromosome}:{start}..{end}:{strand}/{ensembl_assembly_id}'
                   '?content-type=application/json').format(
                chromosome=rnacentral_entry.chromosome.replace('chr', ''),
                start=exon['primary_start'],
                end=exon['primary_end'],
                strand=1 if 'complement' not in exon else -1,
                lncipedia_assembly_id=lncipedia_assembly_id,
                ensembl_assembly_id=ensembl_assembly_id
            )
            ensembl_mapping = requests.get(url)
            if ensembl_mapping.status_code == 200:
                data = ensembl_mapping.json()
                if 'mappings' not in data or len(data['mappings']) == 0:
                    continue
                try:
                    exon['primary_start'] = data['mappings'][0]['mapped']['start']
                    exon['primary_end'] = data['mappings'][0]['mapped']['end']
                    rnacentral_entry.seq_version = '2' # GRCh38
                except: # pylint: disable=W0702
                    print url
                    print data
                    print "Error while mapping to GRCh38"
            else:
                print url
                print 'Ensembl status code: %i' % ensembl_mapping.status_code
            time.sleep(0.1)
        return rnacentral_entry.assembly_info

    def parse_input_data(self):
        """
        Parse the input json file into a Python object.
        """
        with self.input().open('r') as f: # pylint: disable=E1101
            self.data = json.load(f)

    def get_feature_start_end(self, assembly_info):
        """
        """
        primary_start = primary_end = []
        for exon in assembly_info:
            primary_start.append(exon['primary_start'])
            primary_end.append(exon['primary_end'])
        return (min(primary_start), max(primary_end))

    def create_rnacentral_entries(self):
        """
        Process json file into RNAcentralEntry objects that can be written
        to the output files using standard methods.
        """
        for i, seq in enumerate(self.data):
            if self.test and i > 100:
                break

            (feature_location_start, feature_location_end) = self.get_feature_start_end(seq['assembly_info'])

            rnacentral_entry = RNAcentralEntry(
                chromosome=seq['chromosome'],
                common_name=seq['common_name'],
                database='LNCIPEDIA',
                division='XXX', # pylint: disable=W0511
                feature_location_end=feature_location_end,
                feature_location_start=feature_location_start,
                feature_type=seq['feature_type'],
                gene=seq['gene_id'],
                lineage='; '.join(seq['lineage']),
                ncbi_tax_id=seq['ncbi_tax_id'],
                ncrna_class=seq['ncrna_class'],
                optional_id=seq['xref'][2],
                parent_accession=seq['primary_accession'].split('.')[0],
                primary_id=seq['xref'][1],
                product=seq['product'],
                project='PRJ_LNCPD',
                seq_version=seq['primary_accession'].split('.')[-1],
                sequence=seq['sequence'].uppercase(),
                species=seq['scientific_name'],
                references=[
                    {
                        'authors': 'Volders PJ, Verheggen K, Menschaert G, Vandepoele K, Martens L, Vandesompele J, Mestdagh P',
                        'location': 'Nucleic Acids Res. 2015 Jan;43(Database issue):D174-80',
                        'title': 'An update on LNCipedia: a database for annotated human lncRNA sequences',
                        'pmid': 25378313,
                        'doi': '10.1093/nar/gkv295',
                    },
                ],
            )
            for exon in seq['assembly_info']:
                rnacentral_entry.assembly_info.append(exon)

            rnacentral_entry.accession = self.get_accession(rnacentral_entry)
            rnacentral_entry.description = self.get_description(rnacentral_entry)
            if not self.no_assembly_mapping:
                rnacentral_entry.assembly_info = self.map_genomic_assembly(rnacentral_entry)
            self.entries.append(rnacentral_entry)
        print 'Done'

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
        class LncipediaFile(luigi.Task):  # pylint: disable=W0232
            """
            Wrap input file in a luigi task.
            """
            def output(self):
                """
                """
                return luigi.LocalTarget(filename)
        return LncipediaFile()

    def run(self):
        """
        Main task code.
        """
        self.parse_input_data()
        self.create_rnacentral_entries()
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
            outputs[folder] = luigi.LocalTarget(os.path.join(subpath, 'lncipedia_' + folder + '.csv'))
        return outputs


# main entry point
if __name__ == '__main__':
    luigi.run()
