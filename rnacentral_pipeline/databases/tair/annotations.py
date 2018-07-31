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
import collections as coll

from rnacentral_pipeline.databases.helpers import publications as pubs

from rnacentral_pipeline.databases.go_annotations import GoTermAnnotation


ECO_MAPPING = {
    'IBA': '',
    'IC': '',
    'IDA': '',
    'IEA': '',
    'IEP': '',
    'IGI': '',
    'IMP': '',
    'IPI': '',
    'ISM': '',
    'ISS': '',
    'NAS': '',
    'ND': '',
    'RCA': '',
    'TAS': '',
}

TAIR_QUERY = """
select
    acc.locus_tag,
    xrefupi || '_' xref.taxid as rna_id
from rnc_accessions acc
join xref
on
    xref.ac = acc.accession
where
    database = 'TAIR'
    and xref.deleted = 'N'
"""


def load_known_tair(dbconf):
    known = coll.defaultdict(set)
    for result in run_query(dbconf, TAIR_QUERY):
        known[result['locus_tag']].add(result['rna_id'])
    return dict(known)


def publications(raw):
    references = []
    for ref in raw['references'].split('|'):
        prefix, ident = ref.split(':')
        if prefix != 'PMID':
            continue
        references.append(pubs.reference(ident))
    return references


def evidence_code(raw):
    return ECO_MAPPING[raw['eco_term']]


def parse(handle, dbconf):
    reader = csv.DictReader(
        handle,
        fieldnames=[
            'locus',
            'tair_accession',
            'object_name',
            'qualifier',
            'go_term',
            'go_term_id',
            'tair_keyword',
            'aspect',
            'goslim_term',
            'eco_term',
            'eco_description',
            'with',
            'references',
            'source',
            'date',
        ],
        delimiter='\t',
    )
    known = load_known_tair(dbconf)
    for row in reader:
        if row['tair_accession'] not in known:
            continue

        yield GoTermAnnotation(
            rna_id=known[row['tair_accession']],
            qualifier=row['qualifier'].replace(' ', '_'),
            term_id=row['go_term_id'],
            extensions=[],
            evidence_code=evidence_code(row),
            assigned_by=row['source'],
            publications=publications(row),
        )
