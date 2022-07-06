# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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
import operator as op
import typing as ty

from rnacentral_pipeline.databases import data

from . import lookup
from . import helpers

def as_expression(mapping):

    pass




def parse(handle, db_url):
    """
    Loop through a csv file to construct entries from the linked data

    We have to get two columns from the CSV, the URS_taxid, and the Gene ID that
    we should already have joined on to get hits from EA in RNAcentral

    That means there are tuples where there would usually be single values and
    to be able to reuse as much as possible, it might look a bit untidy.

    Also, I build a massive disctionary here, which may not be optimal, for
    looking up the URS -> Gene ID mapping

    """
    rows = csv.DictReader(handle, delimiter=",", quoting=csv.QUOTE_NONE)
    urs_taxids_genes = sorted(map(op.itemgetter("urs_taxid", "GeneID"), rows))

    mappings = lookup.mapping(db_url, urs_taxids_genes)
    gene_lookup = {a[0] : a[1] for a in urs_taxids_genes}
    for urs in urs_taxids_genes:

        yield helpers.as_entry(urs[0], mappings[urs[0]], gene_lookup[urs[0]])
