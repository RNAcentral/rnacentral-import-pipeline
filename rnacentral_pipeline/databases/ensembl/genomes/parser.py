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

import itertools as it
import operator as op
import typing as ty

from Bio import SeqIO

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.ensembl import helpers as common
from rnacentral_pipeline.databases.ensembl.data import Pseudogene
from rnacentral_pipeline.databases.ensembl.genomes import helpers
from rnacentral_pipeline.databases.ensembl.genomes.data import Context
from rnacentral_pipeline.databases.ensembl.vertebrates import helpers as ensembl
from rnacentral_pipeline.databases.helpers import embl


def ncrnas(context: Context, handle) -> ty.Iterable[data.Entry]:
    for record in SeqIO.parse(handle, "embl"):
        current_gene = None
        for feature in record.features:
            if feature.type == "source":
                continue

            if embl.is_gene(feature):
                current_gene = feature
                continue

            if helpers.is_pseudogene(current_gene, feature):
                continue

            if not helpers.is_ncrna(feature):
                continue

            yield helpers.as_entry(context, record, current_gene, feature)


def parse(context: Context, handle) -> ty.Iterable[data.Entry]:
    data = ncrnas(context, handle)
    grouped = it.groupby(data, op.attrgetter("gene"))
    for gene, related in grouped:
        yield from ensembl.generate_related(related)


def pseudogenes(handle: ty.IO) -> ty.Iterable[Pseudogene]:
    try:
        for record in SeqIO.parse(handle, "embl"):
            current_gene = None
            for feature in record.features:
                if feature.type == "source":
                    continue

                if embl.is_gene(feature) and help:
                    current_gene = feature

                if helpers.is_pseudogene(current_gene, feature):
                    gene = embl.gene(feature)
                    if not gene:
                        continue
                    yield Pseudogene(
                        gene=embl.gene(feature),
                        region=common.regions(record, feature)[0],
                    )
    except UnicodeDecodeError:
        import os

        print(f"UTF-8 error in file {handle.name}. Abort parsing.")
        print(f"The working directory is {os.getcwd()}")
        message = f"UFT-8 error in file {handle.name} during pseudogenes parsing. Aborting parse\n"
        message += f"Working directory: {os.getcwd()}"
        slack.send_notification("Ensembl parser error", message)
