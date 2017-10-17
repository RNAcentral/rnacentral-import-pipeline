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

import re

from rnacentral_entry import RNAcentralEntry
from json_parser import JsonParser

from databases.rfam import utils


class RfamSequenceFile(JsonParser):  # pylint: disable=W0232,R0904
    """
    Luigi Task for converting Rfam JSON files into csv files
    that can be loaded into the RNAcentral database.
    """
    database = 'rfam'

    def get_accession(self, entry, database):
        """
        Generate unique accession.
        Example: Z83731.1:29638..30638:rfam
        """
        return ('{parent_accession}.{seq_version}:{start}..{stop}:rfam').format(
            parent_accession=entry.parent_accession,
            seq_version=entry.seq_version,
            start=entry.feature_location_start if not entry.complement else entry.feature_location_end,
            stop=entry.feature_location_end if not entry.complement else entry.feature_location_start,
            external_id=entry.primary_id
        )

    def create_rnacentral_entries(self):
        """
        Process json file into RNAcentralEntry objects that can be written
        to the output files using standard methods.
        """
        skipped = 0
        self.entries = []

        mapping = utils.id_to_insdc_type()
        for i, seq in enumerate(self.data): # pylint: disable=E1101
            if self.test and i > 100: # pylint: disable=E1101
                break

            # replace trailing commas
            seq['lineage'] = re.sub(r'\.$', '', seq['lineage'])

            seq['feature_location_start'] = int(seq['feature_location_start'])
            seq['feature_location_end'] = int(seq['feature_location_end'])

            rnacentral_entry = RNAcentralEntry(
                common_name=seq['common_name'],
                database=self.database.upper(),
                feature_location_end=max(seq['feature_location_start'], seq['feature_location_end']),
                feature_location_start=min(seq['feature_location_start'], seq['feature_location_end']),
                complement=True if seq['feature_location_end'] < seq['feature_location_start'] else False,
                experiment=' '.join(seq['references']),
                feature_type=seq['feature_type'],
                is_composite='N',
                lineage='; '.join(seq['lineage'].split(' ') + [seq['species']]),
                mol_type='full' if seq['is_seed'] == '0' else 'seed',
                ncbi_tax_id=seq['ncbi_tax_id'],
                ncrna_class=mapping.get(seq['primary_id'], None) or seq['ncrna_class'],
                note=' '.join(seq['ontology']) + ' ' + ('Alignment:full' if seq['is_seed'] == '0' else 'Alignment:seed'),
                optional_id=seq['optional_id'],
                parent_accession=seq['parent_accession'],
                primary_id=seq['primary_id'],
                product=seq['description'],
                project='RFAM',
                sequence=seq['sequence'].upper(),
                seq_version=seq['seq_version'],
                species=seq['species'],
                references=[
                    {
                        'authors': 'Nawrocki E.P., Burge S.W., Bateman A., Daub J., Eberhardt R.Y., Eddy S.R., Floden E.W., Gardner P.P., Jones T.A., Tate J., Finn R.D.',
                        'location': 'Nucleic Acids Res. 2015 Jan;43(Database issue):D130-7',
                        'title': 'Rfam 12.0: updates to the RNA families database',
                        'pmid': 25392425,
                        'doi': '10.1093/nar/gku1063',
                    }
                ],
            )

            rnacentral_entry.accession = self.get_accession(rnacentral_entry, self.database) # pylint: disable=E1101
            rnacentral_entry.description = self.get_description(rnacentral_entry) # pylint: disable=E1101
            rnacentral_entry.assembly_info = [{
                'primary_start': rnacentral_entry.feature_location_start,
                'primary_end': rnacentral_entry.feature_location_end,
            }]
            if rnacentral_entry.complement:
                rnacentral_entry.assembly_info[0]['complement'] = ''

            if not rnacentral_entry.is_valid(verbose=True):
                skipped += 1
                continue

            self.entries.append(rnacentral_entry) # pylint: disable=E1101

        print 'Skipped %i sequences' % skipped
