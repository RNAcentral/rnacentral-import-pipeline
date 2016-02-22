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

python path/to/this/file.py
    --local-scheduler
    --destination /path/to/output/files
    --input-file /path/to/input/file.json

Optional parameters:
    --test # to process just the first few entries (default: False)
    --no_assembly_mapping # do not map LNCipedia data onto the latest assembly (default: False)
"""

import luigi
import requests
import time

from rnacentral_entry import RNAcentralEntry
from json_parser import JsonParser


class JsonParserLncipedia(JsonParser):  # pylint: disable=W0232
    """
    Luigi Task for converting LNCipedia Json file into csv files
    that can be loaded into the RNAcentral database.
    """
    no_assembly_mapping = luigi.BoolParameter(default=False, significant=False)
    database = 'lncipedia'

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

    def create_rnacentral_entries(self):
        """
        Process json file into RNAcentralEntry objects that can be written
        to the output files using standard methods.
        """
        for i, seq in enumerate(self.data): # pylint: disable=E1101
            if self.test and i > 100: # pylint: disable=E1101
                break

            (feature_location_start, feature_location_end) = self.get_feature_start_end(seq['assembly_info']) # pylint: disable=E1101

            rnacentral_entry = RNAcentralEntry(
                chromosome=seq['chromosome'],
                common_name=seq['common_name'],
                database=self.database.upper(),
                division='XXX', # pylint: disable=W0511
                feature_location_end=feature_location_end,
                feature_location_start=feature_location_start,
                feature_type=seq['feature_type'],
                gene=seq['gene_id'],
                lineage='; '.join(seq['lineage'] + [seq['scientific_name']]),
                ncbi_tax_id=seq['ncbi_tax_id'],
                ncrna_class=seq['ncrna_class'],
                note=' '.join(seq['ontology']),
                optional_id=seq['xref'][2],
                parent_accession=seq['primary_accession'].split('.')[0],
                primary_id=seq['xref'][1],
                product=seq['product'],
                project='PRJ_LNCPD',
                seq_version=seq['primary_accession'].split('.')[-1],
                sequence=seq['sequence'].upper(),
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

            rnacentral_entry.accession = self.get_accession(rnacentral_entry, self.database) # pylint: disable=E1101
            rnacentral_entry.description = self.get_description(rnacentral_entry) # pylint: disable=E1101
            if not self.no_assembly_mapping:
                rnacentral_entry.assembly_info = self.map_genomic_assembly(rnacentral_entry)
            self.entries.append(rnacentral_entry) # pylint: disable=E1101


# main entry point
if __name__ == '__main__':
    luigi.run(main_task_cls=JsonParserLncipedia)
