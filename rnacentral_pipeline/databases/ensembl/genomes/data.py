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

from pathlib import Path

import attr
from attr.validators import instance_of as is_a

from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.ensembl.gff import load_coordinates


@attr.s()
class Context:
    database = attr.ib(validator=is_a(str), converter=str)
    references = attr.ib(validator=is_a(list))
    gff = attr.ib(validator=is_a(SqliteDict))

    @classmethod
    def build(cls, database: str, references, gff_file: Path) -> "Context":
        return cls(
            database=database,
            references=references,
            gff=load_coordinates(gff_file)
        )

    def accession(self, primary_id: str) -> str:
        return '%s:%s' % (self.database, primary_id)
