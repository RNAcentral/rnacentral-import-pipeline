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

import csv
import sqlite3
import tempfile

import six
import attr
import pytest
import shutil

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.data import IdReference
from rnacentral_pipeline.databases.data import KnownServices
import rnacentral_pipeline.databases.helpers.publications as pub


@pytest.mark.parametrize(
    "raw_id,namespace,external_id",
    [
        (26184978, KnownServices.pmid, "26184978"),
        ("pmid:26184978", KnownServices.pmid, "26184978"),
        ("PMID:26184978", KnownServices.pmid, "26184978"),
        ("doi:10.1038/srep12276", KnownServices.doi, "10.1038/srep12276"),
        ("DOI:10.1038/srep12276", KnownServices.doi, "10.1038/srep12276"),
        ("PMCID:PMC5785218", KnownServices.pmcid, "PMC5785218"),
        ("PMC5785218", KnownServices.pmcid, "PMC5785218"),
        ("pmcid:pmc5785218", KnownServices.pmcid, "PMC5785218"),
        ("pmc5785218", KnownServices.pmcid, "PMC5785218"),
    ],
)
def test_reference_builds_a_id_ref(raw_id, namespace, external_id):
    assert pub.reference(raw_id) == IdReference(
        namespace=namespace,
        external_id=external_id,
    )
