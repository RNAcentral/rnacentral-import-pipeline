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

import csv
import os


class CsvOutputWriter(object): # pylint: disable=W0232
    """
    An object for writing RNAcentralEntry objects into csv files.
    """
    def format_sequence_data(self, entries, fn_short, fn_long):
        """
        Write out the data that is going to be loaded
        into the `load_rnacentral_all` table.
        Short and long sequences need to be written to separate files because
        long sequences are stored as CLOBs in the Oracle database.
        """
        with open(fn_short, 'w') as f_short, \
             open(fn_long, 'w') as f_long:
            writer_short = csv.writer(f_short, delimiter=',', lineterminator='\n')
            writer_long = csv.writer(f_long, delimiter=',', lineterminator='\n')
            for entry in entries:
                row = entry.format_sequence_line()
                if entry.sequence > 4000:
                    writer_long.writerow(row)
                else:
                    writer_short.writerow(row)
        # remove empty files
        if os.stat(fn_short).st_size == 0:
            os.remove(fn_short)
        if os.stat(fn_long).st_size == 0:
            os.remove(fn_long)

    def format_references(self, entries, filename):
        """
        Write out the data for loading into the `load_rnc_references`
        table.
        """
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL, lineterminator='\n')
            for entry in entries:
                writer.writerows(entry.format_references())

    def format_accession_info(self, entries, filename):
        """
        Write out the data for loading into the `load_rnc_accessions` table.
        """
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL, lineterminator='\n')
            for entry in entries:
                writer.writerow(entry.format_ac_line())

    def format_genomic_locations(self, entries, filename):
        """
        Write out the data for loading into the `load_rnc_coordinates` table.
        """
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL, lineterminator='\n')
            for entry in entries:
                writer.writerows(entry.format_genomic_locations())
