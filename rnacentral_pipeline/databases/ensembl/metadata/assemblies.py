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
import itertools as it

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from . import databases as db

REJECTED_TAXIDS = {
    4932,
    5127,
}

BLAT_GENOMES = {
    'anopheles_gambiae',
    'arabidopsis_thaliana',
    'bombyx_mori',
    'caenorhabditis_elegans',
    'dictyostelium_discoideum',
    'drosophila_melanogaster',
    'homo_sapiens',
    'mus_musculus',
    'plasmodium_falciparum',
    'rattus_norvegicus'
    'saccharomyces_cerevisiae',
    'schizosaccharomyces_pombe',
}


EXAMPLE_LOCATIONS = {
    'homo_sapiens': {
        'chromosome': 'X',
        'start': 73819307,
        'end': 73856333,
    },
    'mus_musculus': {
        'chromosome': 1,
        'start': 86351908,
        'end': 86352200,
    },
    'danio_rerio': {
        'chromosome': 9,
        'start': 7633910,
        'end': 7634210,
    },
    'bos_taurus': {
        'chromosome': 15,
        'start': 82197673,
        'end': 82197837,
    },
    'rattus_norvegicus': {
        'chromosome': 'X',
        'start': 118277628,
        'end': 118277850,
    },
    'felis_catus': {
        'chromosome': 'X',
        'start': 18058223,
        'end': 18058546,
    },
    'macaca_mulatta': {
        'chromosome': 1,
        'start': 146238837,
        'end': 146238946,
    },
    'pan_troglodytes': {
        'chromosome': 11,
        'start': 78369004,
        'end': 78369219,
    },
    'canis_familiaris': {
        'chromosome': 19,
        'start': 22006909,
        'end': 22007119,
    },
    'gallus_gallus': {
        'chromosome': 9,
        'start': 15676031,
        'end': 15676160,
    },
    'xenopus_tropicalis': {
        'chromosome': 'NC_006839',
        'start': 11649,
        'end': 11717,
    },
    'saccharomyces_cerevisiae': {
        'chromosome': 'XII',
        'start': 856709,
        'end': 856919,
    },
    'schizosaccharomyces_pombe': {
        'chromosome': 'I',
        'start': 540951,
        'end': 544327,
    },
    'caenorhabditis_elegans': {
        'chromosome': 'III',
        'start': 11467363,
        'end': 11467705,
    },
    'drosophila_melanogaster': {
        'chromosome': '3R',
        'start': 7474331,
        'end': 7475217,
    },
    'bombyx_mori': {
        'chromosome': 'scaf16',
        'start': 6180018,
        'end': 6180422,
    },
    'anopheles_gambiae': {
        'chromosome': '2R',
        'start': 34644956,
        'end': 34645131,
    },
    'dictyostelium_discoideum': {
        'chromosome': 2,
        'start': 7874546,
        'end': 7876498,
    },
    'plasmodium_falciparum': {
        'chromosome': 13,
        'start': 2796339,
        'end': 2798488,
    },
    'arabidopsis_thaliana': {
        'chromosome': 2,
        'start': 18819643,
        'end': 18822629,
    }
}


def reconcile_taxids(taxid):
    """
    Sometimes Ensembl taxid and ENA/Expert Database taxid do not match,
    so to reconcile the differences, taxids in the ensembl_assembly table
    are overriden to match other data.
    """

    taxid = str(taxid)
    if taxid == '284812':  # Ensembl assembly for Schizosaccharomyces pombe
        return 4896  # Pombase and ENA xrefs for Schizosaccharomyces pombe
    return int(taxid)


def is_ignored_assembly(info):
    if info.taxid in REJECTED_TAXIDS:
        return True
    if info.taxid == 7227 and info.division != 'Ensembl':
        return True
    if info.taxid == 5762 and info.assembly_id == 'V1.0':
        return True
    if info.taxid == 6239 and info.division != 'Ensembl':
        return True
    if info.taxid == 8090 and info.ensembl_url != 'oryzias_latipes':
        return True
    return False


class InvalidDomain(Exception):
    """
    Raised when we cannot compute a URL for a domain.
    """
    pass


@attr.s()
class AssemblyExample(object):
    chromosome = attr.ib(validator=is_a(basestring))
    start = attr.ib(validator=is_a(int))
    end = attr.ib(validator=is_a(int))

    @classmethod
    def build(cls, raw, example_locations):
        key = raw['species.url'].lower()
        example = example_locations.get(key, None)
        if not example:
            return None

        return cls(
            chromosome=str(example['chromosome']),
            start=example['start'],
            end=example['end'],
        )


@attr.s()
class AssemblyInfo(object):
    assembly_id = attr.ib(validator=is_a(basestring))
    assembly_full_name = attr.ib(validator=is_a(basestring))
    gca_accession = attr.ib(validator=optional(is_a(basestring)))
    assembly_ucsc = attr.ib(validator=optional(is_a(basestring)))
    common_name = attr.ib(validator=optional(is_a(basestring)))
    taxid = attr.ib(validator=is_a(int))
    ensembl_url = attr.ib(validator=is_a(basestring))
    division = attr.ib(validator=is_a(basestring))
    blat_mapping = attr.ib(validator=is_a(bool))
    example = attr.ib(validator=optional(is_a(AssemblyExample)))

    @classmethod
    def build(cls, raw, example_locations):
        url = raw['species.url'].lower()
        is_mapped = url in BLAT_GENOMES

        return cls(
            assembly_id=raw['assembly.default'],
            assembly_full_name=raw['assembly.name'],
            gca_accession=raw.get('assembly.accession', None),
            assembly_ucsc=raw.get('assembly.ucsc_alias', None),
            common_name=raw.get('species.common_name', None),
            taxid=reconcile_taxids(raw['species.taxonomy_id']),
            ensembl_url=url,
            division=raw['species.division'],
            blat_mapping=is_mapped,
            example=AssemblyExample.build(raw, example_locations),
        )

    @property
    def subdomain(self):
        """Given E! division, returns E!/E! Genomes url."""

        if self.division == 'Ensembl':
            return 'ensembl.org'
        if self.division == 'EnsemblPlants':
            return 'plants.ensembl.org'
        if self.division == 'EnsemblMetazoa':
            return 'metazoa.ensembl.org'
        if self.division == 'EnsemblBacteria':
            return 'bacteria.ensembl.org'
        if self.division == 'EnsemblFungi':
            return 'fungi.ensembl.org'
        if self.division == 'EnsemblProtists':
            return 'protists.ensembl.org'
        if self.division == 'EnsemblVertebrates':
            return 'vertebrates.ensembl.org'
        raise InvalidDomain(self.division)

    def writeable(self):
        chromosome = None
        start = None
        end = None
        if self.example:
            chromosome = self.example.chromosome
            start = self.example.start
            end = self.example.end

        return [
            self.assembly_id,
            self.assembly_full_name,
            self.gca_accession,
            self.assembly_ucsc,
            self.common_name,
            self.taxid,
            self.ensembl_url,
            self.division,
            self.subdomain,
            chromosome,
            start,
            end,
            int(self.blat_mapping),
        ]


def fetch(connections, query_handle, example_locations):
    results = db.run_queries_across_databases(connections, query_handle)
    for (_, rows) in results:
        raw = {r['meta_key']: r['meta_value'] for r in rows}
        if raw['species.division'] == 'EnsemblBacteria':
            continue
        info = AssemblyInfo.build(raw, example_locations)
        if is_ignored_assembly(info):
            continue
        yield info


def write(connections, query, output):
    """
    Parse the given input handle and write the readable data to the CSV.
    """

    data = fetch(connections, query, EXAMPLE_LOCATIONS)
    data = it.imap(op.methodcaller('writeable'), data)
    csv.writer(output).writerows(data)
