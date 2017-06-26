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

from pgloader import PGLoader
from rfam.hits_csv import RfamHitsCSV
from rfam.pgload_families import PGLoadFamilies
from rfam.pgload_clans import PGLoadClans
from rfam.config import db as DBConfig

CONTROL_FILE = """
LOAD CSV
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
  INTO postgresql://{user}:{password}@{host}:{port}/{db}?load_rnc_rfam_model_hits
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
     WITH skip header = 0,
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
    (load.sequence_stop - load.sequence_start)::float / rna.len,
    load.model_start,
    load.model_stop,
    (load.model_stop - load.model_start)::float / models.length,
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
"""


class PGLoadHits(PGLoader):
    def requires(self):
        return [
            RfamHitsCSV(),
            PGLoadFamilies(),
            PGLoadClans(),
        ]

    def filename(self):
        return self.requires()[0].output().fn

    def control_file(self):
        """
        This generates the text to write as a control file for pgloader.

        :returns: The control file text
        """
        config = DBConfig()
        return CONTROL_FILE.format(
            filename=self.filename,
            user=config.user,
            password=config.password,
            host=config.host,
            port=config.port,
            db=config.db_name,
        )


if __name__ == '__main__':
    luigi.run(main_task_cls=PGLoadHits)
