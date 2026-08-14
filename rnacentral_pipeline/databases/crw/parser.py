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

Metadata comes from R2DT's crw-metadata.tsv rather than r2dt_models: that table
carries neither taxid nor cellular_location -- files/r2dt/load-models.ctl drops
both on the way in -- so the query this used to run could never have worked.
Reading the repo also means the import no longer depends on the analyze stage
having populated a table first.
"""

import typing as ty
from pathlib import Path

from Bio import SeqIO

from rnacentral_pipeline.databases.crw import helpers
from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.rnacentral.r2dt.models import crw as models


def parse(metadata_handle: ty.IO, sequences: Path) -> ty.Iterable[Entry]:
    indexed = SeqIO.index(str(sequences), "fasta")
    for info in models.parse(metadata_handle):
        entry = helpers.as_entry(info, indexed)
        if entry:
            yield entry
