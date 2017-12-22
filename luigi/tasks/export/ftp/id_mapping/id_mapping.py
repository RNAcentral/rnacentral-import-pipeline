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
import csv

import luigi

from tasks.config import export


SQL = """
SELECT
    xref.upi,
    xref.ac AS accession,
    xref.taxid,
    acc.external_id,
    acc.optional_id,
    acc.feature_name,
    acc.ncrna_class,
    acc.gene,
    db.descr AS database
FROM xref,
     rnc_accessions acc,
     rnc_database db
WHERE
    rna.upi = xref.upi
    AND xref.ac = acc.accession
    AND xref.dbid = db.id
    AND xref.deleted = 'N'
"""


class IdMapping(luigi.Task):
    def output(self):
        return luigi.LocalTarget(export().id_mapping('id_mapping.tsv'))

    def mappings(self, sql=SQL):
        with cursor() as cur:
            cur.execute(sql)
            for result in cur:
                database = result['database']
                if database == 'PDBE':
                    database = 'PDB'
                gene = result['gene'] or ''
                gene = gene.replace('\t', ' ')
                accession = result['external_id']
                if database == 'ENA' or database == 'HGNC':
                    accession = result['accession']

                yield [
                    result['upi'],
                    database,
                    accession,
                    result['taxid'],
                    result['rna_type'],
                    gene,
                ]

    def run(self):
        try:
            os.makedirs(self.output().fn)
        except:
            pass

        with open(self.output().fn, 'w') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerows(self.mappings())
