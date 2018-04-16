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

import csv

import attr
from attr.validators import instance_of as is_a

from rnacentral.psql import PsqlWrapper

PRECUSOR_MAPPING_QUERY = """
SELECT
    precursor_xref.upi as precursor,
    mature_xref.upi as mature,
    precursor_xref.taxid
FROM xref precursor_xref
JOIN rnc_accessions precursor_acc
ON
    precursor_xref.ac = precursor_acc.accession
    AND precursor_acc.feature_name = 'precursor_RNA'
JOIN rnc_accessions mature_acc
ON
    mature_acc.external_id = precursor_acc.external_id
    AND mature_acc.feature_name != 'precursor_RNA'
JOIN xref mature_xref
ON
    mature_xref.ac = mature_acc.accession
    AND mature_xref.dbid = precursor_xref.dbid
    AND precursor_xref.upi != mature_xref.upi
    AND precursor_xref.taxid = mature_xref.taxid
    AND mature_xref.ac != precursor_acc.accession
WHERE
    precursor_xref.dbid = 4
    AND precursor_xref.deleted = 'N'
    AND mature_xref.deleted = 'N'
"""


QUERY = """
WITH mirna_precursors as (
{precursor_query}
)
SELECT
    pre.upi,
    pre.taxid,
    max(pre.description),
    max(pre.rna_type),
    array_agg(mirna_precursors.precursor) as precursor_rnas
FROM rnc_rna_precomputed pre
LEFT JOIN mirna_precursors ON mirna_precursors.mature = pre.upi
WHERE
    pre.taxid IS NOT NULL
    AND rna_type IS NOT NULL
    AND description IS NOT NULL
group by pre.upi, pre.taxid
""".format(
    precursor_query=PRECUSOR_MAPPING_QUERY,
)

FIELDS = [
    'database',
    'DB_Object_ID',
    'DB_Object_Symbol',
    'DB_Object_Name',
    'DB_Object_Synonym',
    'DB_Object_Type',
    'Taxon',
    'Parent_Object_ID',
    'DB_Xref',
    'Gene_Product_Properties',
]


class GpiEntry(object):
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=is_a(int))
    description = attr.ib(validator=is_a(basestring))
    so_term = attr.ib(validator=is_a(basestring))
    precursor_rnas = attr.ib(validator=is_a(list), default=[])

    @classmethod
    def build(cls, result):
        return cls(
            upi=result['upi'],
            taxid=result['taxid'],
            description=result['description'],
            so_term=result['so_term'],
        )

    @property
    def rna_id(self):
        return '{upi}_{taxid}'.format(upi=self.upi, taxid=self.taxid)

    def precursor_info(self):
        return 'precursor_rna={upis}'.format(
            upis=','.join(self.precursor_rnas),
        )

    def as_dict(self):
        return {
            'database': 'RNAcentral',
            'DB_Object_ID': self.rna_id,
            'DB_Object_Symbol': '',
            'DB_Object_Name': self.description,
            'DB_Object_Synonym': '',
            'DB_Object_Type': self.so_term,
            'Taxon': 'taxon:%s' % self.taxid,
            'Parent_Object_ID': '',
            'DB_Xref': '',
            'Gene_Product_Properties': self.precursor_info(),
        }


def entries(db):
    psql = PsqlWrapper(db)
    for result in psql.copy_to_iterable(QUERY):
        entry = GpiEntry.build(result)
        yield entry.as_dict()


def export(db, handle):
    handle.write('!gpi-version: 1.2\n')
    writer = csv.DictWriter(handle, FIELDS)
    writer.writerows(entries(db))
