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

import os
import json

from jsonschema import validate

from rnacentral.psql import PsqlWrapper


MOD_URL = 'http://modomics.genesilico.pl/sequences/list/{id}'

SCHEMA = os.path.realpath(os.path.join(
    os.path.dirname(__file__),
    '..',
    '..',
    '..',
    '..',
    'files',
    'ensembl-schema.json'
))

BASE_SQL = """
SELECT
    json_build_object(
        'rnacentral_id', pre.id,
        'description', max(pre.description),
        'sequence', max(case
            when rna.seq_short is null then rna.seq_long
            else rna.seq_short
        end),
        'md5', max(rna.md5),
        'rna_type', max(pre.rna_type),
        'taxon_id', max(xref.taxid),
        'xrefs', array_agg(
            json_build_object(
                'database', db.display_name,
                'external_id', acc.external_id,
                'optional_id', acc.optional_id,
                'molecule_type', acc.mol_type
            )
        )
    )
FROM rna
JOIN xref ON xref.upi = rna.upi
JOIN rnc_rna_precomputed pre ON pre.upi = xref.upi AND pre.taxid = xref.taxid
JOIN rnc_database db ON db.id = xref.dbid
JOIN rnc_accessions acc
ON
    xref.ac = acc.accession
WHERE
    xref.deleted = 'N'
    AND %s
group by pre.id
"""

SINGLE_SQL = BASE_SQL % "xref.upi = '{upi}' and xref.taxid = {taxid}"
RANGE_SQL = BASE_SQL % "rna.id BETWEEN {min_id} AND {max_id}"

TRUSTED_DB = set([
    'gtrnadb',
    'lncrnadb',
    'mirbase',
    'modomics',
    'pdbe',
    'snopy'
    'srpdb',
    'tmrna website',
])


def external_id(data):
    if data['database'] == 'PDBe':
        return '%s_%s' % (data['external_id'], data['optional_id'])
    if data['database'] == 'Modomics':
        return MOD_URL.format(id=data['external_id'])
    return data['external_id']

def is_high_quality(data):
    name = data['database'].lower()
    if name in TRUSTED_DB:
        return True
    if name == 'rfam':
        return data['molecule_type'] == 'seed'
    return False


def as_xref(xref):
    return {'database': xref['database'], 'id': external_id(xref)}


def builder(data):
    result = dict(data)
    xrefs = []
    for xref in data['xrefs']:
        if is_high_quality(xref):
            xrefs.append(as_xref(xref))

    result['xrefs'] = xrefs
    result['sequence'] = result['sequence'].replace('T', 'U')
    return result


def export(db, query, **kwargs):
    psql = PsqlWrapper(db)
    for result in psql.copy_to_iterable(query, **kwargs):
        try:
            data = json.loads(result['json_build_object'])
            yield builder(data)
        except:
            raise


def upi(db, upi, taxid):
    results = export(db, SINGLE_SQL, upi=upi, taxid=taxid)
    try:
        return next(results)
    except StopIteration:
        raise ValueError("Found no entries for %s_%i" % (upi, taxid))


def range_filename(min_id, max_id):
    return 'ensembl-xrefs-{min}-{max}.json'.format(min=min_id, max=max_id)


def range(db, min_id, max_id):
    return list(export(db, RANGE_SQL, min_id=min_id, max_id=max_id))


def write_and_validate(handle, results, schema_file=SCHEMA):
    with open(schema_file, 'r') as raw:
        schema = json.load(raw)

    validate(results, schema)
    json.dump(handle, results)
