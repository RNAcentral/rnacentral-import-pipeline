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

from rnacentral.psql import PsqlWrapper


QUERY = """
select
    pre.id,
    go_terms.go_term_id,
    string_agg(hits.rfam_model_id, '|')
from rnc_rna_precomputed pre
join xref on xref.upi = pre.upi and xref.taxid = pre.taxid
join rfam_model_hits hits on hits.upi = pre.upi
join rfam_go_terms go_terms on hits.rfam_model_id = go_terms.rfam_model_id
where
    pre.taxid is not null
    and xref.deleted = 'N'
    and pre.rfam_problems is not null
    and pre.rfam_problems != ''
    and '{"has_issue": false}'::jsonb <@ pre.rfam_problems::jsonb
group by pre.id, go_terms.go_term_id
"""


def export(config, handle):
    psql = PsqlWrapper(config)
    psql.write_query(handle, QUERY, use='tsv')
