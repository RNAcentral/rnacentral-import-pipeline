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

import luigi

from tasks.config import db
from .utils.db import cursor


TABLES = [
    'load_rnacentral',
    'load_rnacentral_all',
    'load_rnc_accessions',
    'load_rnc_references',
    'load_rnc_coordinates',
    'load_upi_max_versions',
    'load_rnc_secondary_structure',
    'load_retro_tmp',
    'load_max_versions',
]


class CleanupRelease(luigi.Task):  # pylint: disable=R0904
    """
    This will empty out the load_* tables. This is distinct from the normal
    loading procedure because sometimes things go wrong it and it is very hard
    to debug when all the intermediate data has been deleted.
    """

    def run(self):
        with cursor(db()) as cur:
            for table in TABLES:
                cur.execute('truncate table %s' % table)
                cur.commit()
