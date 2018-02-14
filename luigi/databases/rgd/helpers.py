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

import re
import csv
import gzip
import tempfile
import operator as op
import itertools as it
import collections as coll
from contextlib import contextmanager

import attr
from Bio import SeqIO

from databases.data import Exon
from databases.data import Entry
from databases.helpers import phylogeny as phy


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
basic_rna_type = op.itemgetter('GENE_TYPE')
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

    @classmethod
    def from_name(cls, name):
        return cls(
            organism=name,
        )

    @property
    def sequence_path(self):
        if self.organism == 'rat':
            return 'pub/data_release/UserReqFiles/RAT_NUCLEOTIDE_SEQUENCES.fa.gz'
        raise ValueError("Cannot find sequences for %s" % self.organism)

    @property
    def gene_path(self):
        filename = 'GENES_%s.txt' % (self.organism.upper())
        return 'pub/data_release/{filename}'.format(filename=filename)

    def sequence_uri(self, config):
        return 'ftp://{host}/{path}'.format(
            host=config.host,
            path=self.sequence_path,
        )

    def gene_uri(self, config):
        return 'ftp://{host}/{path}'.format(
            host=config.host,
            path=self.gene_path,
        )


def corrected_records(handle):
    """
    This just strips out the ',' that is at the end of ids that RGD provides.
    """

    seen = coll.defaultdict(set)
    for record in SeqIO.parse(handle, 'fasta'):

        if not str(record.seq):
            continue

        # These are probably protein, so skip them
        if record.id.startswith('XM_') or record.id.startswith('NM_'):
            continue

        # Change given ids into a probably unique id
        given = record.id.replace(',', '')
        match = re.search('gene RGD:(\d+),', record.description)
        if not match:
            raise ValueError("RGD fasta must state gene id")

        record.id = '{given}-{gene}'.format(
            given=given,
            gene=match.group(1),
        )

        # Prevent writing duplicate entries
        if str(record.seq) in seen[record.id]:
            continue

        seen[record.id].add(str(record.seq))
        yield record


@contextmanager
def indexed(filename):
    """
    This will decompress the fasta file we fetch from RGD and process it to
    produce an indexed file that uses the gene id as the key to index the
    sequences.
    """

    with tempfile.NamedTemporaryFile() as tmp:
        with gzip.open(filename, 'r') as raw:
            SeqIO.write(corrected_records(raw), tmp, 'fasta')

        tmp.flush()
        yield SeqIO.index(tmp.name, 'fasta')


def rna_type(entry):
    basic = basic_rna_type(entry)
    if basic == 'rrna':
        return 'rRNA'
    if basic == 'trna':
        return 'tRNA'
    if basic == 'ncrna':
        return 'other'
    if basic == 'snrna':
        return 'snRNA'
    raise ValueError("Cannot determine rna type for %s" % entry)


def indexed_gene_id(entry):
    return primary_id(entry)


def accession(entry, index=None):
    if index is not None:
        return 'RRID:RGD_%s:%i' % (primary_id(entry), index)
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
    return basic_rna_type(entry) in KNOWN_RNA_TYPES


def url(entry):
    return URL.format(id=primary_id(entry))


def seq_version(_):
    return '1'


def seq_xref_ids(entry):
    """
    This will produce the list of all ids that could be used to extract the
    """

    xref_ids = []
    for ids in xref_data(entry).values():
        xref_ids.extend(ids)
    return xref_ids


def sequences_for(entry, sequences):
    count = 0
    seqs = set()
    for xref_id in seq_xref_ids(entry):
        record_id = xref_id + '-' + primary_id(entry)
        if record_id not in sequences:
            continue
        count += 1
        record = sequences[record_id]
        seqs.add(str(record.seq))

    # if not seqs:
    #     raise ValueError("No sequences found for: %s" % entry)

    # assert count == len(seqs), "RGD is a pita"
    return list(seqs)


def exons(entry):
    if not entry['CHROMOSOME_6.0']:
        return []
    complement = False
    strand = entry['STRAND_6.0']
    if strand == '-':
        complement = True
    return [Exon(
        chromosome='chr%s' % entry['CHROMOSOME_6.0'],
        primary_start=int(entry['START_POS_6.0']),
        primary_end=int(entry['STOP_POS_6.0']),
        complement=complement,
    )]


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
    if entry['OLD_SYMBOL']:
        return entry['OLD_SYMBOL'].split(';')
    return []


def as_entry(data, seqs):
    sequences = sequences_for(data, seqs)
    for index, sequence in enumerate(sequences):
        acc_index = index + 1
        if len(sequences) == 1:
            acc_index = None

        return Entry(
            primary_id=primary_id(data),
            accession=accession(data, acc_index),
            ncbi_tax_id=taxid(data),
            database='RGD',
            sequence=sequence,
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
