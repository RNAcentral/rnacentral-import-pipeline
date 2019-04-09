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

from .data import Interactor, Interaction


def identifiers(raw):
    return [Identifier.build(r) for r in raw.split('|')]


def interactor(row, offset):
    parts = {
        'unique_ids': row[offset],
        'alternatives': row[offset + 2],
        'aliases': row[offset + 4],
        'taxids': row[offset + 9],
        'xrefs': row[offset + 22]
    }

    unique_ids = identifiers(parts['unique_ids'])
    assert len(unique_ids) == 1

    taxids = identifiers(parts['taxid'])
    assert len(taxids) > 0
    taxid = taxids[0].value

    return Interactor(
        identifiers=unique_ids[0],
        taxonomy_id=tax_id,
        alternatives=identifiers(parts['identifiers']),
        aliases=identifiers(parts['aliases']),
        xrefs=identifiers(parts['xrefs']),
    )


def interactor_a(row):
    return interactor(row, 0)


def interactor_b(row):
    return interactor(row, 1)


def involves_rnacentral(row):
    return interactor_a(row).database == 'rnacentral' or \
            interactor_b(row).database == 'rnacentral'


def interactors(row):
    int_a = interactor_a(row)
    int_b = interactor_b(row)
    if a.database == 'rnacentral':
        return int_a, int_b
    if b.database == 'rnacentral':
        return int_b, int_a
    raise ValueError("Could not determine which involves RNacentral")


def interactions(row):
    return None


def publications(row):
    pmids = six.moves.map(lambda p: p.key == 'pubmed', identifiers(row[8]))
    return [reference(pmid.value) for pmid in pmids]


def as_annotation(row):
    urs, interactor = interactors(row)
    pubs = publications(row)
    for interaction in interactions(row):
        yield Interaction(
            interactor=interactor,
            urs=urs,
            interaction=interaction,
            publications=pubs,
        )
