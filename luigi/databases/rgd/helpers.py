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

import csv
import operator as op
import itertools as it

import attr

from databases.data import Entry
from databases import helpers as phy


KNOWN_RNA_TYPES = set([
    'ncrna',
    'rrna',
    'snrna',
    'trna',
])

URL = 'https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={id}'

XREF_NAMING = {
    'ensembl': 'ENSEMBL_ID',
    'ncbi_gene': 'NCBI_GENE_ID',
    'genbank': 'GENBANK_NUCLEOTIDE',
}

primary_id = op.itemgetter('GENE_RGD_ID')
rna_type = op.itemgetter('GENE_TYPE')
pmids = op.itemgetter('CURATED_REF_RGD_ID', 'CURATED_REF_PUBMED_ID',
                      'UNCURATED_PUBMED_ID')
gene = op.itemgetter('SYMBOL')
locus_tag = op.itemgetter('SYMBOL')


def known_organisms():
    """
    Get the names of organisms that RGD annotates. This will only return 'rat'
    currently. While they have annotations for other databases we only process
    rat because this is the only organism that they provide chromosomes that we
    can extract sequences from. Without that I would have to write logic, like
    in Rfam, that tracks down genomes. This isn't needed yet.
    """
    return [
        'rat'
    ]


@attr.s()
class RgdInfo(object):
    organism = attr.ib()
    genome = attr.ib()

    @classmethod
    def from_name(cls, name):
        genome = None
        if name == 'rat':
            genome = 'rn6'
        elif name == 'human':
            return 'hg38'
        else:
            raise ValueError("Cannot determine genome for %s" % name)

        return cls(
            organism=name,
            genome=genome,
        )

    @property
    def pretty(self):
        return self.organism[0].upper() + self.organism[1:]

    @property
    def chromosome_path(self):
        return 'pub/data_release/agr/fasta/{genome}'.format(genome=self.genome)

    @property
    def gene_path(self):
        filename = 'GENES_%s.txt' % (self.organism.upper())
        return 'pub/data_release/{filename}'.format(filename=filename)

    @property
    def gff_path(self):
        return 'data_release/GFF3/Gene/{name}'.format(name=self.pretty)

    def chromosomes(self, conn):
        return conn.lst(self.chromosome_path())

    def gff_files(self, conn):
        paths = conn.lst(self.gff_path())
        return [path for path in paths if 'RATMINE' not in path]


def accession(entry):
    return 'RRID:RGD_%s' % primary_id(entry)


def taxid(_):
    return 10116


def species(entry):
    return phy.species(taxid(entry))


def fetch_and_split(entry, name):
    if not entry[name]:
        return None
    return entry[name].split(';')


def is_ncrna(entry):
    return rna_type(entry) in KNOWN_RNA_TYPES


def url(entry):
    return URL.format(id=primary_id(entry))


def seq_version(_):
    return '1'


def sequence(entry, seqs):
    return ''
    # record = seqs[primary_id(entry)]
    # return str(record.seq)


def exons(_):
    return []


def xref_data(entry):
    xrefs = {}
    for name, key in XREF_NAMING.items():
        data = fetch_and_split(entry, key)
        if data:
            xrefs[name] = data
    return xrefs


def references(entry):
    refs = []
    # refs.append(pub.as_reference(25355511))  # The general RGD citation
    possible_ids = it.chain.from_iterable(pmids(entry))
    for idset in possible_ids:
        pubids = idset.split(';')
        for _ in pubids:
            # refernces.append(pub.as_reference
            pass
    return refs


def description(entry):
    tag = None
    if locus_tag(entry):
        tag = ' (%s)' % locus_tag(entry)

    return '{species} {name}{tag}'.format(
        name=entry['NAME'],
        species=species(entry),
        tag=tag,
    )


def gene_synonyms(entry):
    if entry['OLD_NAME']:
        return entry['OLD_SYMBOL'].split(';')
    return []


def as_entry(data, seqs):
    return Entry(
        primary_id=primary_id(data),
        accession=accession(data),
        ncbi_tax_id=taxid(data),
        database='RGD',
        sequence=sequence(data, seqs),
        exons=exons(data),
        rna_type=rna_type(data),
        url=url(data),
        seq_version=seq_version(data),

        xref_data=xref_data(data),

        gene=gene(data),
        locus_tag=locus_tag(data),
        gene_synonyms=gene_synonyms(data),
        description=description(data),

        references=references(data),
    )


def as_rows(lines):
    return csv.DictReader(lines, delimiter='\t')
