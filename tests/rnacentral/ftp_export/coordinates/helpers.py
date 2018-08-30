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

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import data

from tests.helpers import run_with_replacements


def fetch_raw(rna_id, assembly):
    path = os.path.join('files', 'ftp-export', 'genome_coordinates',
                        'query.sql')
    _, taxid = rna_id.split('_')
    return run_with_replacements(
        path,
        (':assembly_id', "'%s'" % assembly),
        (':taxid', taxid),
        ('WHERE\n',
         '''where
         pre.id = '%s'
         and''' % rna_id
         ), take_all=True)


def fetch_coord(rna_id, assembly):
    return data.parse(fetch_raw(rna_id, assembly))


def fetch_all(taxid, assembly):
    path = os.path.join('files', 'ftp-export', 'genome_coordinates',
                        'query.sql')
    return data.parse(run_with_replacements(
        path,
        (':assembly_id', "'%s'" % assembly),
        (':taxid', str(taxid)),
        take_all=True,
    ))
