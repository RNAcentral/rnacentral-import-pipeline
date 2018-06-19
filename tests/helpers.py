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
import tempfile

from rnacentral_pipeline.psql import PsqlWrapper


def run_range_as_single(upi, path):
    parts = upi.split('_')
    upi, taxid = parts[0], int(parts[1])
    with tempfile.NamedTemporaryFile('w') as tmp:
        # path = os.path.join('files', 'search-export', 'query.sql')
        with open(path, 'rb') as raw:
            query = raw.read()
            query = query.replace(
                'rna.id BETWEEN :min AND :max',
                "xref.upi ='{upi}' AND xref.taxid = {taxid}".format(
                    upi=upi,
                    taxid=taxid,
                )
            )
            tmp.write(query)
            tmp.flush()
        psql = PsqlWrapper(os.environ['PGDATABASE'])
        results = psql.copy_file_to_iterable(tmp.name, json=True)

        try:
            return next(results)
        except StopIteration:
            raise ValueError("Found no entries for %s_%i" % (upi, taxid))
