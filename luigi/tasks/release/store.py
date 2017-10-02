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
from .prepare import PrepareRelease
from .utils.db import cursor


SQL = """
select
    dbid,
    id,
    release_type,
    release_date,
    force_load
from rnacen.rnc_release
where status = 'L'
order by id
"""

CREATE_INDEX = """
create index if not exists load_rnacentral_all$database
on rnacen.load_rnacentral_all(database)
"""


class StoreRelease(luigi.Task):  # pylint: disable=R0904
    """
    This will store all the loaded data. It will not prepare a new database
    release by itself.
    """

    def run(self):
        with cursor(db()) as cur:
            cur.exectue(CREATE_INDEX)
            cur.execute("select rnc_update.prepare_releases('F')")
            cur.execute(SQL)
            for result in cur.fetchall():
                cur.execute('select rnc_update.new_update_release(%s, %s)',
                            (result[0], result[1]))
                cur.commit()
