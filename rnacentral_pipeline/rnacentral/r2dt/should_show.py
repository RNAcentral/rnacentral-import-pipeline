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
import operator as op

from rnacentral_pipeline.rnacentral.r2dt import data


def parse(handle: ty.IO) -> ty.Iterable[data.ShowInfo]:
    for record in csv.DictReader(handle):
        yield data.ShowInfo.from_raw(record)


def write(handle: ty.IO, output: ty.IO):
    writer = csv.writer(output)
    data = parse(handle)
    rows = map(op.methodcaller("writeable"), data)
    writer.writerows(rows)
