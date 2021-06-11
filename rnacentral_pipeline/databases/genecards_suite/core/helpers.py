# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.helpers import phylogeny as phy

from . import data


def primary_id(context: data.Context, row) -> str:
    gene = context.gene(row)
    rna_id = context.urs(row)
    return "%s:%s:%s" % (context.database, gene, rna_id)


def accession(context: data.Context, row) -> str:
    return primary_id(context, row)


def taxid(context, row) -> int:
    rna_id = context.urs(row)
    return int(rna_id.split("_", 1)[1])


def species(context: data.Context, row) -> str:
    return phy.species(taxid(context, row))


def lineage(context: data.Context, row) -> str:
    return phy.lineage(taxid(context, row))


def common_name(context: data.Context, row) -> str:
    return phy.common_name(taxid(context, row))
