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


QUERY = """
select
    hits.upi,
    hits.rfam_model_id,
    score,
    e_value,
    sequence_start,
    sequence_stop,
    model_start,
    model_stop,
    models.long_name
from rfam_model_hits hits
join rna_active active on active.upi = hits.upi
join rfam_models models on models.rfam_model_id = hits.rfam_model_id
order by hits.upi, hits.sequence_start, hits.rfam_model_id
"""


def write(connection, out):
    command = "COPY ({query}) to STDOUT".format(
        query=QUERY.replace('\n', ' '),
    )
    cursor = connection.cursor()
    cursor.copy_expert(command, out)
