# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import typing as ty
import csv

from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source

RFAM_QUERY = """
select 
	rfam_model_id,
	'rfam',
	131567,
	so_rna_type,
	rna_type 
from rfam_models 
where
	so_rna_type is not null
"""

def load_info(db_url: str) -> ty.Dict[str, ty.Tuple[str, int]]:
    return {}


def parse(cm_stat: ty.IO, db_url: str) -> ty.Iterable[ModelInfo]:
    known_info = load_info(db_url)
    for row in csv.reader(cm_stat):
        info = known_info[row[0]]
        yield ModelInfo(
            model_name=row[0],
            so_rna_type=info[0],
            taxid=info[1],
            source=Source.rfam,
            length=int(row[1]),
            basepairs=int(row[2]),
        )
