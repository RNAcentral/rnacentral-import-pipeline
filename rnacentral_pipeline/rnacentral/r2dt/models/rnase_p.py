# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import csv
import typing as ty

from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source


def parse(handle, extra=None) -> ty.Iterable[ModelInfo]:
    for row in csv.DictReader(handle, delimiter="\t"):
        so_term_id = "SO:0000386"
        taxid = int(row["taxid"])
        yield ModelInfo(
            model_name=row["model_name"],
            so_rna_type=so_term_id,
            taxid=taxid,
            source=Source.rnase_p,
            cell_location=None,
            length=None,
            basepairs=None,
        )
