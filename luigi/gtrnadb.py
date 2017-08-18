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

import json
from hashlib import md5

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from databases import helpers as dbs

from ensembl.data import Exon


def optionally(instance_type, **kwargs):
    """
    Return an attribute that is either none or of the given type.
    """
    return attr.ib(
        validator=optional(is_a(instance_type)),
        default=attr.Factory(None),
        **kwargs
    )


def possibly_empty(instance_type, **kwargs):
    """
    Return an attribute that defaults to being empty and must be of the given
    type.
    """
    return attr.ib(
        validator=is_a(instance_type),
        default=attr.Factory(instance_type),
        **kwargs
    )


def url(data):
    """
    Adds a http:// to the start of the url in data.
    """
    return 'http://%s' % data['url']


def anticodon(data):
    """
    Get the anticodon of this entry.
    """
    return data['metadata']['anticodon']


def note_data(entry):
    """
    Create the dict that will be stored as a note. This is basically whatever
    GtRNAdb gives us, with a few duplicate or unneeded fields removed.
    """

    note = {}
    note.update(entry)
    del note['organism']
    del note['pseudogene']
    return note


def chromosome(location):
    """
    Get the chromosome this location is part of.
    """

    chrom = location['chromosome']
    if chrom == 'Chromosome':
        return '1'
    return chrom


def common_name(data):
    """
    Get a standardized common name for the given taxon id.
    """
    return dbs.common_name(data['ncbi_tax_id'])


def lineage(data):
    """
    Get a standardized lineage for the given taxon id.
    """
    return dbs.lineage(data['ncbi_tax_id'])


def species(data):
    """
    Get a standardized species name for the given taxon id.
    """
    return dbs.species(data['ncbi_tax_id'])


def description(data):
    """
    Generate a description for the entries specified by the data.
    """
    return '{name} tRNA-{aa} ({anticodon})'.format(
        name=common_name(data),
        aa=data['metadata']['isotype'],
        anticodon=anticodon(data),
    )


def product(data):
    """
    Generate the product for the entries specified by the data.
    """
    return 'tRNA-{aa} {anticodon}'.format(
        aa=data['metadata']['isotype'],
        anticodon=data['metadata']['anticodon'],
    )


def primary_id(data, location):
    """
    Generate a primary key for the given data and location.
    """

    start = min(e['start'] for e in location['exons'])
    stop = max(e['stop'] for e in location['exons'])
    return '{gene}:{accession}:{start}-{stop}'.format(
        gene=data['gene'],
        accession=location['exons'][0]['INSDC_accession'],
        start=start,
        stop=stop,
    )


@attr.s(frozen=True)
class SecondaryStructure(object):
    """
    This represents the secondary structure from GtRNAdb.
    """
    dot_bracket = attr.ib(validator=is_a(basestring))

    @classmethod
    def empty(cls):
        """
        Create an empty secondary structure.
        """
        return cls(dot_bracket='')

    @classmethod
    def from_gtrnadb(cls, raw):
        """
        Generate a secondary structure from the raw angle bracket string. This
        will transform it into a reasonable dot-bracket string and create a
        SecondaryStructure object.
        """

        transformed = raw.\
            replace('>', '(').\
            replace('<', ')')
        assert set(transformed) == set('(.)')
        return cls(dot_bracket=transformed)

    def __bool__(self):
        """
        Check if this is empty.
        """
        return bool(self.dot_bracket)

    @property
    def md5(self):
        """
        Compute the MD5 of the dot_bracket string.
        """
        return md5(self.dot_bracket).hexdigest()


def exons_from_gtrnadb(data):
    """
    This will create the Exons from the data provided by GtRNAdb.
    """

    exons = []
    for exon in data['locations']:
        complement = None
        if exon['strand'] == '+':
            complement = False
        elif exon['strand'] == '-':
            complement = True
        else:
            raise ValueError("Invalid strand %s" % exon)

        exons.append(Exon(
            primary_start=int(exon['start']),
            primary_end=int(exon['stop']),
            complement=complement,
        ))
    return exons


@attr.s(frozen=True)
class Entry(object):
    """
    This represents an RNAcentral entry from GtRNAdb that we will write out for
    import.
    """

    primary_id = attr.ib(validator=is_a(basestring))
    accession = attr.ib(validator=is_a(basestring))
    ncbi_tax_id = attr.ib(validator=is_a(int))
    database = attr.ib(validator=is_a(basestring))
    sequence = attr.ib(validator=is_a(basestring))
    exons = attr.ib(validator=is_a(list))
    rna_type = attr.ib(validator=is_a(basestring))
    url = attr.ib(validator=is_a(basestring))

    note_data = possibly_empty(dict)
    xref_data = possibly_empty(dict)

    chromosome = optionally(str)
    species = optionally(str)
    common_name = optionally(str)
    lineage = optionally(str)
    gene = optionally(str)
    locus_tag = optionally(str)
    optional_id = optionally(str)
    product = optionally(str)
    parent_accession = optionally(str)
    ordinal = optionally(str)
    non_coding_id = optionally(str)
    project = optionally(str)
    keywords = optionally(str)
    division = optionally(str)
    organelle = optionally(str)
    allele = optionally(str)
    anticodon = optionally(str)
    experiment = optionally(str)
    function = optionally(str)
    inference = optionally(str)
    map = optionally(str)
    old_locus_tag = optionally(str)
    operon = optionally(str)
    standard_name = optionally(str)
    description = optionally(str)

    gene_synonyms = possibly_empty(list)
    references = possibly_empty(list)

    secondary_structure = attr.ib(
        default=attr.Factory(SecondaryStructure.empty),
        validator=is_a(SecondaryStructure),
    )

    @property
    def db_xrefs(self):
        """
        Return a JSON encoded dict representing the xref data.
        """
        return json.dumps(self.xref_data)

    @property
    def note(self):
        """
        Return a JSON encoded dictionary representing the note data.
        """
        return json.dumps(self.note_data)

    @property
    def feature_type(self):
        """
        Return the feature for the RNA type.
        """
        if self.rna_type in set(['rRNA', 'tRNA', 'precursor_RNA', 'tmRNA']):
            return 'misc_RNA'
        return 'ncRNA'

    @property
    def ncrna_class(self):
        """
        The ncRNA class. If the feature type is not ncRNA this this will be the
        empty string.
        """
        if self.feature_type != 'ncRNA':
            return ''
        return self.rna_type

    def gene_synonym(self):
        """
        Returns a comma separated list of gene synonyms.
        """
        return ','.join(self.gene_synonyms)

    @classmethod
    def from_gtrnadb(cls, data):
        """
        Take an entry from GtRNAdb and produce the RNAcentrals that it
        represents. A single entry may represent more than one Entry because it
        may occur in more than one location. As we provide an accession for
        each location this ends up representing more than one RNAcentral Entry.
        """

        if data['metadata']['pseudogene']:
            return

        two_d = SecondaryStructure.from_gtrnadb(data['secondary_structure'])
        for location in data['genome_locations']:
            parent_ac = location['exons'][0]['INSDC_accession']
            accession = '{ac}:{gene}'.format(ac=parent_ac, gene=data['gene'])

            yield cls(
                primary_id=primary_id(data, location),
                accession=accession,
                ncbi_tax_id=int(data['ncbi_tax_id']),
                database='GtRNADB',
                sequence=data['sequence'],
                exons=exons_from_gtrnadb(location),
                rna_type='tRNA',
                url=url(data),
                note_data=note_data(data),
                secondary_structure=two_d,
                chromosome=chromosome(location),
                species=species(location),
                common_name=common_name(data),
                anticodon=anticodon(data),
                lineage=lineage(location),
                gene=data['gene'],
                optional_id=data['gene'],
                product=product(data),
                parent_accession=parent_ac,
                description=description(data),
                mol_type='genomic DNA',
                feature_location_start=1,
                feature_location_stop=len(data['sequence']),
                gene_synonyms=data.get('synonyms', []),
            )


def parse(filename):
    """
    This will parse a JSON file produced by GtRNAdb and yield the RNAcentral
    entries that it represents.
    """

    with open(filename, 'rb') as raw:
        data = json.load(raw)
        for entry in data:
            for result in Entry.from_gtrnadb(entry):
                yield result
