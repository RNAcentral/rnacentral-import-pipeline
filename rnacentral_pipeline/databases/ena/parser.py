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

import typing as ty
import logging
from pathlib import Path

from Bio import SeqIO

import rnacentral_pipeline.databases.helpers.embl as embl
from rnacentral_pipeline.databases.data import Entry

from rnacentral_pipeline.databases.ena import context, dr, helpers, ribovore
from rnacentral_pipeline.databases.ena import mapping as tpa

LOGGER = logging.getLogger(__name__)


class InvalidEnaFile(Exception):
    """
    This is raised when there is something wrong with the ENA EMBL file.
    """
    pass


def parse(ctx: context.Context, path: Path) -> ty.Iterable[Entry]:
    """
    Parse a file like object into an iterable of Entry objects. This will parse
    each feature in all records of the given EMBL formatted file to produce the
    Entry objects.
    """

    for record in SeqIO.parse(str(path), "embl"):
        if len(record.features) == 0:
            LOGGER.warn("Skipping record %s with no features" % record.id)
            continue

        ctx.add_total()
        if len(record.features) != 2:
            raise InvalidEnaFile(
                "ENA EMBL files must have 2 features/record %s" % record
            )

        feature = record.features[1]
        if helpers.is_protein(feature):
            LOGGER.info("Skipping mis-annotated protein: %s", record.id)
            ctx.add_skipped_protein()
            continue

        if helpers.is_pseudogene(feature):
            LOGGER.info("Skipping pseudogene")
            ctx.add_skipped_pseudogene()
            continue

        if record.id not in ctx.dr:
            raise InvalidEnaFile("Somehow parsed DR refs are for wrong record")

        entry = helpers.as_entry(ctx, record, feature)
        if ctx.ribovore is not None:
            ribo_result = ctx.ribovore.get(record.id, None)
            if helpers.is_skippable_sequence(entry, ribo_result):
                LOGGER.info(
                    f"Skipping record ({record.id}) excluded by ribotyper")
                ctx.add_riboytper_skip()
                continue

        ctx.add_parsed()
        yield entry


def parse_with_context(ctx: context.Context, path: Path) -> ty.Iterable[Entry]:
    entries = parse(ctx, path)
    return ctx.expand_tpa(entries)
