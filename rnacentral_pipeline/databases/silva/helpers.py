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
"""

import re
import operator as op

from rnacentral_pipeline.databases import data
import rnacentral_pipeline.databases.helpers.publications as pubs
import rnacentral_pipeline.databases.helpers.phylogeny as phy

url = op.itemgetter('silvaUri')

CLASS_PATTERN = re.compile(r';\s*$')

KNOWN_TYPES = {
    'rRNA': 'SO:0000252',
    'rRNA_12S': 'SO:0002128',
    'rRNA_16S': 'SO:0001000',
    'rRNA_18S': 'SO:0000407',
    'small_subunit_rRNA': 'SO:0000650',
}


def primary_id(row):
    return 'SILVA:%s:%s' % (row['insdcAccession'], row['location'])


def taxid(row):
    return int(row['ncbiTaxId'])


def sequence(row):
    return row['sequence'].replace('U', 'T')


def inference(row):
    value = row['classification']
    value = value.replace(';', '; ')
    return re.sub(CLASS_PATTERN, '', value)


def version(row):
    _, version = row['insdcAccession'].split('.', 1)
    return version


def rna_type(row):
    given = row['type']
    if given in KNOWN_TYPES:
        return KNOWN_TYPES[given]
    raise ValueError("Unknown RNA type")


def lineage(row):
    return phy.lineage(taxid(row))


def species(row):
    return phy.species(taxid(row))


def common_name(row):
    return phy.common_name(taxid(row))


def as_entry(row):
    try:
        return data.Entry(
            primary_id=primary_id(row),
            accession=primary_id(row),
            ncbi_tax_id=taxid(row),
            database='SILVA',
            sequence=sequence(row),
            regions=[],
            rna_type=rna_type(row),
            url=url(row),
            seq_version=version(row),
            common_name=common_name(row),
            species=species(row),
            lineage=lineage(row),
            references=[
                pubs.reference('doi:10.1093/nar/gks1219'),
            ],
            inference=inference(row),
        )
    except phy.UnknownTaxonId:
        return None
