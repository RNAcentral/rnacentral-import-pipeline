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

RiboVision used to be imported by scraping rnacentral_mapping.html from
apollo.chemistry.gatech.edu, which keyed entries on PDB chains. That host is
gone and RiboVision2 serves a JavaScript app with no equivalent page, so the
models are taken from the copy R2DT ships instead, the same way CRW works.

The metadata lives in the repo rather than r2dt_models because that table
carries neither taxid nor cellular_location -- files/r2dt/load-models.ctl drops
both on the way in.
"""

import typing as ty
from pathlib import Path

from Bio import SeqIO

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.ribovision import helpers
from rnacentral_pipeline.rnacentral.r2dt.models import ribovision as models


def parse(directory: Path, sequences: Path) -> ty.Iterable[Entry]:
    """
    Read every metadata.tsv under an R2DT data directory. RiboVision splits its
    models across a large and a small subunit directory.
    """
    indexed = SeqIO.index(str(sequences), "fasta")
    for metadata in sorted(directory.glob("ribovision-*/metadata.tsv")):
        with metadata.open("r") as handle:
            for info in models.parse(handle):
                entry = helpers.as_entry(info, indexed)
                if entry:
                    yield entry
