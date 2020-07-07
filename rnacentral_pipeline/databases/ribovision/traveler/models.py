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

import csv
import operator as op

from rnacentral_pipeline.rnacentral.traveler.data import Source
from rnacentral_pipeline.rnacentral.traveler.data import ModelInfo


def lookup_taxid(species):
    if species == 'Deinococcus radiodurans':
        return 1299
    if species == 'Bacillus subtilis':
        return 1423
    if species == 'Spinacia oleracea':
        return 3562
    if species == 'Dictyostelium discoideum':
        return 44689
    if species == 'Drosophilla melanogaster':
        return 7227
    if species == 'Escherichia coli':
        return 562
    if species == 'Haloarcula marismortui':
        return 2238
    if species == 'Homo sapiens':
        return 9606
    if species == 'Kluyveromyces lactis':
        return 28985
    if species == 'Mycolicibacterium smegmatis':
        return 1772
    if species == 'Mycobacterium tuberculosis':
        return 1773
    if species == 'Plasmodium falciparum':
        return 5833
    if species == 'Staphylococcus aureus':
        return 1280
    if species == 'Saccharomyces cerevisiae':
        return 4932
    if species == 'Tetrahymena thermophila':
        return 5911
    if species == 'Thermus thermophilus':
        return 274
    if species == 'Drosophila melanogaster':
        return 7227
    if species == 'Trypanosoma brucei':
        return 5691
    raise ValueError("Unknown species name: " + species)


def as_location(raw):
    if not raw:
        return None
    if raw == 'mito':
        return 'Mitochondrion'
    if raw == 'chloroplast':
        return 'Chloroplast'
    raise ValueError("Unknown raw location: " + raw)


def so_term(row):
    if row['model_name'] == 'DR_LSU_3D':
        return 'SO:0001001'
    if 'LSU' in row['model_name'] or '23S' in row['model_name'] or \
            '28S' in row['model_name'] or row['model_name'] == 'F3H4_G_18380' \
            or row['model_name'] == 'GC14_75':
        return 'SO:0000651'
    if 'SSU' in row['model_name']:
        return 'SO:0000650'
    raise ValueError("Could not figure out SO term for: %s" % row)


def parse(handle):
    for row in csv.DictReader(handle, delimiter='\t'):
        so_term_id = so_term(row)

        if not row['taxid']:
            taxid = lookup_taxid(row['species'])
        else:
            taxid = int(row['taxid'])

        location = as_location(row['cellular_location'])
        yield ModelInfo(
            model_id=row['model_name'],
            is_intronic=False,
            so_term=so_term_id,
            taxid=taxid,
            accessions=[],
            source=Source.ribovision,
            cell_location=location,
            length=None,
        )


def write(handle, output):
    data = parse(handle)
    data = map(op.methodcaller('writeable'), data)
    csv.writer(output).writerows(data)
