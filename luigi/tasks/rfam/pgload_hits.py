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

from tasks.utils.pgloader import PGLoader

from tasks.rfam.hits_csv import RfamHitsCSV
from tasks.rfam.pgload_families import RfamPGLoadFamilies

CONTROL_FILE = """LOAD CSV
     FROM '{filename}'
     HAVING FIELDS
      (
      UPI,
      SEQUENCE_START,
      SEQUENCE_STOP,
      STRAND,
      RFAM_MODEL_ID,
      MODEL_START,
      MODEL_STOP,
      OVERLAP,
      E_VALUE,
      SCORE
      )
  INTO {db_url}?load_rnc_rfam_model_hits
     TARGET COLUMNS
      (
      UPI,
      SEQUENCE_START,
      SEQUENCE_STOP,
      RFAM_MODEL_ID,
      MODEL_START,
      MODEL_STOP,
      OVERLAP,
      E_VALUE,
      SCORE
      )
     WITH skip header = 1,
          fields escaped by double-quote,
          fields terminated by ','
BEFORE LOAD DO
$$
truncate table load_rnc_rfam_model_hits;
$$
AFTER LOAD DO
$$ insert into rnc_rfam_model_hits (
    sequence_start,
    sequence_stop,
    sequence_completness,
    model_start,
    model_stop,
    model_completeness,
    overlap,
    e_value,
    score,
    rfam_model_id,
    upi
) (
select
    load.sequence_start,
    load.sequence_stop,
    (load.sequence_stop - load.sequence_start)::float / rna.len::float,
    load.model_start,
    load.model_stop,
    (load.model_stop - load.model_start)::float / models.length::float,
    load.overlap,
    load.e_value,
    load.score,
    load.rfam_model_id,
    load.upi
from rna, rnc_rfam_models as models, load_rnc_rfam_model_hits as load
where
    rna.upi = load.upi
    and models.rfam_model_id = load.rfam_model_id
);
$$,
$$
truncate table load_rnc_rfam_model_hits;
$$
;
"""


class RfamPGLoadHits(PGLoader):  # pylint: disable=R0904
    """
    This task will run pgloader on the generate hits CSV file. This will
    populate the database with the hits. The loading process does not attempt
    to prevent duplicates so reimporting the same data will lead to duplicates
    in the database.
    """

    def requires(self):
        return [
            RfamHitsCSV(),
            RfamPGLoadFamilies(),
        ]

    def control_file(self):
        filename = RfamHitsCSV().output().fn
        return CONTROL_FILE.format(filename=filename)
