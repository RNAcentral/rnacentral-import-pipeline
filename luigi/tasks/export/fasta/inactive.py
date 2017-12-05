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

from tasks.config import export

from .utils import FastaExportBase


CREATE_INACTIVE_TABLE = """
CREATE TEMP TABLE rna_inactive as
select upi from rna
"""

DELETE_FROM_INACTIVE_TABLE = """
delete from rna_inactive
using rna_active as active
where
    active.upi = rna_inactive.upi;
"""

FETCH_INACTIVE = """
select
    pre.id,
    pre.description,
    case when rna.seq_short is null
        then rna.seq_long
        else rna.seq_short end as sequence
from rna_inactive inactive
join rnc_rna_precomputed pre on pre.upi = inactive.upi
join rna on rna.upi = inactive.upi
where
    pre.id = rna.upi
    and pre.upi = inactive.upi
"""


class InactiveFastaExport(FastaExportBase):
    table = CREATE_INACTIVE_TABLE
    populate = DELETE_FROM_INACTIVE_TABLE
    fetch = FETCH_INACTIVE

    def output(self):
        return luigi.LocalTarget(export().ftp(
            'sequences',
            'rnacentral_inactive.fasta',
        ))
