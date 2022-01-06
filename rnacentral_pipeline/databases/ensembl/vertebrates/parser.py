# -*- coding: utf-8 -*-

"""
Copyright [2010-2018] EMBL-European Bioinformatics Institute
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

import itertools as it
import logging
import operator as op
import typing as ty
from pathlib import Path

import attr
from Bio import SeqIO

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import embl
from rnacentral_pipeline.databases.ensembl.data import TranscriptInfo
from rnacentral_pipeline.databases.ensembl.gencode import helpers as gencode

from rnacentral_pipeline.databases.ensembl import helpers as common
from rnacentral_pipeline.databases.ensembl.data import Pseudogene
from rnacentral_pipeline.databases.ensembl.vertebrates import helpers
from rnacentral_pipeline.databases.ensembl.vertebrates.context import Context

LOGGER = logging.getLogger(__name__)

IGNORE_FEATURES = {
    "source",
    "STS",
    "misc_feature",
}


def as_entry(
    record, gene, feature, context: Context, is_nonchromosomal=False
) -> ty.Optional[data.Entry]:
    """
    Turn the Record, Gene feature, transcript feature and Context into a Entry
    object for output.
    """

    species, common_name = helpers.organism_naming(record)
    xref_data = embl.xref_data(feature)

    try:
        sequence = embl.sequence(record, feature)
    except Exception as err:
        LOGGER.exception(err)
        return None

    pid = helpers.primary_id(feature)
    if pid in context.gff:
        info = context.gff[pid]
    else:
        if is_nonchromosomal:
            info = TranscriptInfo.from_feature(record, feature)
        else:
            LOGGER.warn("Cannot find transcript info for %s: %s", pid, feature)
            return None

    entry = data.Entry(
        primary_id=pid,
        accession=helpers.accession(feature),
        ncbi_tax_id=embl.taxid(record),
        database="ENSEMBL",
        sequence=sequence,
        regions=info.regions,
        rna_type=info.so_rna_type,
        url=helpers.url(feature),
        seq_version=helpers.seq_version(feature),
        lineage=embl.lineage(record),
        chromosome=helpers.chromosome(record),
        parent_accession=record.id,
        common_name=common_name,
        species=species,
        gene=embl.gene(gene),
        locus_tag=embl.locus_tag(gene),
        optional_id=embl.gene(gene),
        note_data=helpers.note_data(feature),
        xref_data=xref_data,
        product=helpers.product(feature),
        references=helpers.references(),
        mol_type="genomic DNA",
        pseudogene="N",
        is_composite="N",
    )

    return attr.evolve(entry, description=helpers.description(context, gene, entry))


def ncrnas(raw, context: Context) -> ty.Iterable[data.Entry]:
    """
    This will parse an EMBL file for all Ensembl Entries to import.
    """

    is_nonchromosomal = "nonchromosomal" in raw.name
    for record in SeqIO.parse(raw, "embl"):
        current_gene = None
        for feature in record.features:

            if feature.type in IGNORE_FEATURES:
                LOGGER.debug("Skipping ignored feature type for %s", feature)
                continue

            if embl.is_gene(feature):
                current_gene = feature
                continue

            if helpers.is_pseudogene(current_gene, feature):
                LOGGER.debug("Skipping psuedogene %s", feature)
                continue

            if not helpers.is_ncrna(feature):
                LOGGER.debug("Skipping feature %s because it is not ncRNA", feature)
                continue

            entry = as_entry(
                record,
                current_gene,
                feature,
                context,
                is_nonchromosomal=is_nonchromosomal,
            )
            if not entry:
                LOGGER.warn("Could not parse %s" % feature)
                continue

            if context.is_supressed(entry):
                LOGGER.debug("Skipping supressed Rfam family %s", feature)
                continue

            if context.is_excluded(entry):
                LOGGER.debug("Skipping excluded entry %s", feature)
                continue

            yield entry


def parse(
    raw: ty.IO, gff_file: Path, family_file=None, excluded_file=None
) -> ty.Iterable[data.Entry]:
    """
    This will parse an EMBL file for all Ensembl Entries to import.
    """

    context = Context.build(
        gff_file, family_file=family_file, excluded_file=excluded_file
    )
    loaded = ncrnas(raw, context)
    grouped = it.groupby(loaded, op.attrgetter("gene"))
    for _, entries in grouped:
        related = list(entries)
        for entry in helpers.generate_related(related):
            yield entry

        from_gencode = []
        for entry in related:
            if context.from_gencode(entry):
                from_gencode.append(gencode.update_entry(entry))
        for entry in helpers.generate_related(from_gencode):
            yield entry


def pseudogenes(handle: ty.IO) -> ty.Iterable[Pseudogene]:
    for record in SeqIO.parse(handle, "embl"):
        current_gene = None
        for feature in record.features:
            if feature.type in IGNORE_FEATURES:
                LOGGER.debug("Skipping ignored feature type for %s", feature)
                continue

            if embl.is_gene(feature):
                current_gene = feature

            if helpers.is_pseudogene(current_gene, feature):
                gene = embl.gene(feature)
                if not gene:
                    continue
                yield Pseudogene(
                    gene=gene,
                    region=common.regions(record, feature)[0],
                )
