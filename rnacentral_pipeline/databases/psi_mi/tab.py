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

import re
import csv
import logging
import typing as ty
from datetime import date
import collections as coll

import tatsu

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import publications as pubs

LOGGER = logging.getLogger(__name__)

IDENT_GRAMMAR = r'''
@@grammar::PSI_MI_IDENTIFIER

start = empty:empty | idents $ ;

empty = '-' ;

idents = '|'.{ ident }+ ;

ident = xref:string ':' [ value:string [ '(' description:string ')' ] ] ;

string = ('"' quoted_string  '"') | simple_string ;

simple_string = /[^|():\t]+/ ;

quoted_string = /(\\"|[^"])+/ ;
'''

class IdentSemantics(object):
    def string(self, ast):
        if isinstance(ast, (list, tuple)):
            l, v, r = ast
            ast = v
        ast = ast.replace(r'\"', '"')
        return ast


IDENT_MODEL = tatsu.compile(IDENT_GRAMMAR, semantics=IdentSemantics())


def identifiers(raw: str) -> ty.List[data.InteractionIdentifier]:
    assert raw, "Must have at least one identifier"
    idents = []
    for possible in IDENT_MODEL.parse(raw):
        if 'empty' in possible:
            continue
        idents.append(data.InteractionIdentifier(
            possible['xref'], 
            possible.get('value', ''),
            possible.get('description', None)
        ))

    return idents


def as_taxid(value):
    if value == '-':
        return None
    taxids = identifiers(value)
    unique = set()
    for tid in taxids:
        if tid.value not in unique:
            unique.add(tid.value)
    assert len(unique) == 1, "Too many taxids in: %s" % value
    return int(taxids[0].value)


def as_unique_id(value):
    unique_ids = identifiers(value)
    if not unique_ids:
        return None
    assert len(unique_ids) <= 1, "Too many ids: %s" % value
    if unique_ids:
        return unique_ids[0]
    return None


def as_bool(value):
    if value.lower() == 'true':
        return True
    if value.lower() == 'false':
        return False
    raise ValueError("Unknown bool value: %s" % value)


def as_date(value):
    year, month, day = value.split('/')
    return date(int(year), int(month), int(day))


def as_pubs(value):
    refs = []
    for ident in identifiers(value):
        if ident.key == 'pubmed':
            refs.append(pubs.reference(ident.value))
    return refs


def stoichiometry(value):
    if value == '-':
        return None
    return int(value)


def as_interactor(row, interactor: data.InteractorType) -> ty.Optional[data.Interactor]:
    field_names = {
        'ID(s) interactor': ('id', as_unique_id),
        'Alt. ID(s) interactor': ('alt_ids', identifiers),
        'Alias(es) interactor': ('aliases', identifiers),
        'Taxid interactor': ('taxid', as_taxid),
        'Biological role(s) interactor': ('biological_role', identifiers),
        'Experimental role(s) interactor': ('experimental_role', identifiers),
        'Type(s) interactor': ('interactor_type', identifiers),
        'Xref(s) interactor': ('xrefs', identifiers),
        'Annotation(s) interactor': ('annotations', identifiers),
        'Feature(s) interactor': ('features', identifiers),
        'Stoichiometry(s) interactor': ('stoichiometry', stoichiometry),
        'Identification method participant': ('participant_identification',
                                              identifiers),
    }

    parts: ty.Dict[str, ty.Any] = {}
    for field_template, (key, fn) in field_names.items():
        if interactor == data.InteractorType.A and \
                field_template == 'ID(s) interactor':
            field_template = '#ID(s) interactor'

        field_name = '%s %s' % (field_template, interactor.name)
        parts[key] = fn(row[field_name])

    if not parts['id']:
        return None
    return data.Interactor(**parts)


def as_interaction(row) -> ty.Optional[data.Interaction]:
    field_names = {
        'Interaction detection method(s)': ('methods', identifiers),
        'Publication Identifier(s)': ('publications', as_pubs),
        'Interaction type(s)': ('types', identifiers),
        'Source database(s)': ('source_database', identifiers),
        'Interaction identifier(s)': ('ids', identifiers),
        'Confidence value(s)': ('confidence', identifiers),
        'Expansion method(s)': ('methods', identifiers),
        'Interaction Xref(s)': ('xrefs', identifiers),
        'Interaction annotation(s)': ('annotations', identifiers),
        'Host organism(s)': ('host_organisms', as_taxid),
        'Creation date': ('create_date', as_date),
        'Update date': ('update_date', as_date),
        'Negative': ('is_negative', as_bool),
    }

    parts: ty.Dict[str, ty.Any] = {}
    for field_name, (key, fn) in field_names.items():
        parts[key] = fn(row[field_name])

    parts['interactor1'] = as_interactor(row, data.InteractorType.A)
    parts['interactor2'] = as_interactor(row, data.InteractorType.B)
    
    if not parts['interactor1'] or not parts['interactor2']:
        return None

    return data.Interaction(**parts)


def parse(handle) -> ty.Iterator[data.Interaction]:
    rows = csv.DictReader(handle, delimiter='\t')
    interactions = map(as_interaction, rows)
    valid = filter(None, interactions)
    return valid
