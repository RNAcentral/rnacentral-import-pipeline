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
import itertools as it

import attr

def empty_to_none(raw):
    if not raw:
        return None
    return raw


@attr.s()
class RfamDatabaseLink(object):
    rfam_family = attr.ib()
    database = attr.ib()
    comment = attr.ib(convert=empty_to_none)
    external_id = attr.ib()
    other = attr.ib(convert=empty_to_none)

    @classmethod
    def from_row(cls, row):
        """
        Build an object from a dictionary that has the same names as the
        database columns.
        """

        external_id = row['db_link']
        if row['db_id'] in {'SO', 'GO'}:
            external_id = '%s:%s' % (row['db_id'], row['db_link'])

        return cls(
            rfam_family=row['rfam_acc'],
            database=row['db_id'],
            comment=row['comment'],
            external_id=external_id,
            other=row['other_params'],
        )


def parse(handle):
    reader = csv.DictReader(handle, delimiter='\t')
    return it.imap(RfamDatabaseLink.from_row, reader)
