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

import csv
import operator as op
import collections as coll

from Bio import SeqIO
from Bio import Entrez

from rnacentral_pipeline.databases.helpers import embl
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pub


ALLOWED_RNA = {
    "miscRNA",
    "ncRNA",
    "rRNA",
    "scRNA",
    "snRNA",
    "snoRNA",
    "tRNA",
}

row_rna_type = op.itemgetter("type_of_gene")


def value(row, name, required=False):
    current = row[name]
    if current == "-":
        current = None
    if required:
        assert current is not None, "Value missing for: %s" % name
    return current


def gene_id(row):
    return value(row, "GeneID", required=True)


def taxid(row):
    return int(value(row, "#tax_id", required=True))


def row_is_ncrna(row):
    return row_rna_type(row) in ALLOWED_RNA


def ncrna_feature(row):
    ncrna = None
    for feature in row["sequence"].features:
        if feature.type in {"ncRNA", "misc_RNA"}:
            if ncrna is not None:
                raise ValueError("Multiple ncRNAs")
            ncrna = feature
    return ncrna


def rna_type(row):
    rna_type = row_rna_type(row)
    feature = ncrna_feature(row)
    if not feature:
        return rna_type
    return embl.rna_type(feature)


def primary_id(row):
    return "NCBI_GENE:" + gene_id(row)


def accession(row):
    return "NCBI_GENE:" + gene_id(row)


def sequence(row):
    return str(row["sequence"].seq)


def seq_version(_):
    return "1"


def url(row):
    return "https://www.ncbi.nlm.nih.gov/gene/" + gene_id(row)


def xref_data(row):
    data = coll.defaultdict(list)
    given = value(row, "dbXrefs") or ""
    if not given:
        return {}
    for dbid in given.split(","):
        parts = dbid.split(":", 1)
        assert len(parts) == 2, "Invalid DB id: " + dbid
        data[db].append(key)
    return data


def gene(row):
    return value(row, "Symbol", required=True)


def locus_tag(row):
    return value(row, "LocusTag")


def gene_synonyms(row):
    raw = value(row, "Synonyms")
    if not raw:
        return []
    return raw.split(",")


def references(row):
    return [pub.reference(25355515)]


def description(row):
    template = [species(row), gene(row)]
    name = common_name(row)
    if name:
        name = "(%s)" % name
        template.insert(1, name)
    return " ".join(template)


def species(row):
    return phy.species(taxid(row))


def lineage(row):
    return phy.lineage(taxid(row))


def common_name(row):
    return phy.common_name(taxid(row))


def product(row):
    feature = ncrna_feature(row)
    if feature:
        return embl.product(feature)
    return None


def ncrnas(handle):
    found = False
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row_is_ncrna(row):
            found = True
            yield row

    if not found:
        raise ValueError("No ncRNAs found")
