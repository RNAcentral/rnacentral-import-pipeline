# -*- coding: utf-8 -*-

"""
Copyright [2009-present] EMBL-European Bioinformatics Institute
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
import logging
import typing as ty

from rnacentral_pipeline.databases.ensembl import urls
from rnacentral_pipeline.databases.ensembl.urls import GenesetUrls

LOGGER = logging.getLogger(__name__)


def rows(
    selected: ty.Iterable[GenesetUrls], assemblies: ty.IO
) -> ty.Iterable[ty.List[str]]:
    by_gca = {s.accession: s for s in selected if s.genome_url}
    skipped = 0
    for species, assembly_id, taxid, gca in csv.reader(assemblies):
        found = by_gca.get(gca)
        if not found:
            skipped += 1
            continue
        yield [species, assembly_id, taxid, found.genome_url, found.gff_url]
    if skipped:
        LOGGER.warning("Skipped %i assemblies with no genome in species.json", skipped)


def write(location: str, assemblies: ty.IO, output: ty.IO):
    csv.writer(output).writerows(rows(urls.urls_for(location), assemblies))
