# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.db import connection


def get_active_sequence_count(cur):
    """
    Uses a similar query to the one that gets active sequences to count the number of
    sequences that should be in the parquet file.
    """
    query = """select
  count(pc.upi)
  from rnc_rna_precomputed pc
  join rna on rna.upi = pc.upi
    where is_active and taxid is not NULL
    """
    cur.execute(query)
    active_seq_count = cur.fetchone()[0]

    return active_seq_count
