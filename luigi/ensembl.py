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
import re
import csv
import logging
import itertools as it

from rnacentral_entry import RNAcentralEntry

import attr
import luigi
from Bio import SeqIO

from attr.validators import instance_of as is_a


logger = logging.getLogger(__name__)


@attr.s()
class Output(object):
    short_sequences = attr.ib(validator=is_a(luigi.LocalTarget))
    long_sequences = attr.ib(validator=is_a(luigi.LocalTarget))
    references = attr.ib(validator=is_a(luigi.LocalTarget))
    accession = attr.ib(validator=is_a(luigi.LocalTarget))
    locations = attr.ib(validator=is_a(luigi.LocalTarget))

    @classmethod
    def target(cls, name):
        pass

    @classmethod
    def build(cls, base, database):
        return cls(
        )


@attr.s()
class OutputHandles(object):
    short_sequences = attr.ib(validator=is_a(file))
    long_sequences = attr.ib(validator=is_a(file))
    references = attr.ib(validator=is_a(file))
    accessions = attr.ib(validator=is_a(file))
    locations = attr.ib(validator=is_a(file))

    @classmethod
    def build(cls, output):

        return cls(
            short_sequences=open(output.short_sequences, 'w'),
            long_sequences=open(output.long_sequences, 'w'),
            references=open(output.references, 'w'),
            accessions=open(output.accesion, 'w'),
            locations=open(output.locations, 'w'),
        )

    def __enter__(self):
        def as_csv(handle):
            return csv.writer(handle, delimiter=',', quotechar='"',
                              quoting=csv.QUOTE_ALL, lineterminator='\n')

        return OutputHandles(
            short_sequences=as_csv(self.short_sequences),
            long_sequences=as_csv(self.long_sequences),
            references=as_csv(self.references),
            accessions=as_csv(self.accessions),
            locations=as_csv(self.locations),
        )

    def __exit__(self, *args):
        self.short_sequences.close()
        self.long_sequences.close()
        self.references.close()
        self.accessions.close()
        self.locations.close()


@attr.s(slots=True, frozen=True)
class GeneInfo(object):
    locus_tag = attr.ib()
    description = attr.ib()

    @classmethod
    def from_feature(cls, feature):
        return cls(
            locus_tag=feature.qualifiers.get('locus_tag', [None])[0],
            description=feature.qualifiers.get('note', [None])[0],
        )


def qualifier_value(feature, name, pattern, max_allowed=1):
    values = set()
    for note in feature.qualifiers.get(name, []):
        match = re.match(pattern, note)
        if match:
            values.add(match.group(1))
    if max_allowed is not None and len(values) > max_allowed:
        raise ValueError("Multiple values (%s) for %s",
                         ', '.join(sorted(values)), name)
    if len(values) == 0:
        return None
    if max_allowed == 1:
        return values.pop()
    return values


class EnsemblImporter(luigi.Task):

    def format(self):
        return 'embl'

    def write(self, writers, data):
        if data.is_short():
            writers.short_sequences.writerow(data.format_sequence())
        else:
            writers.long_sequences.writerow(data.format_sequence())

        writers.references.writerow(data.format_references())
        writers.accessions.writerow(data.format_accession())
        writers.locations.writerow(data.format_locations())

    def output(self):
        if not os.path.exists(self.destination):
            os.makedirs(self.destination)
        return Output.build(self.destination, 'ensembl')

    def gene_info(self, records):
        info = {}
        for record in records:
            for feature in record.features:
                if feature.type != 'gene':
                    continue
                info[feature.id] = GeneInfo.from_feature(feature)
        return info

    def data(self):
        data = SeqIO.parse(self.filename(), self.format())
        data = it.chain.from_iterable(data.features)
        # Maybe should process all genes to find the names of things? That way
        # we can have a better description than something super generic based
        # off the annotations? We could also extract the locus tags this way.
        data = it.ifilterfalse(data, self.is_pseduogene)
        data = it.ifilter(data, self.is_ncrna)
        data = it.imap(self.create_rnacentral_entries, data)
        data = it.ifilter(data, lambda e: e.is_valid())

    def run(self):
        with OutputHandles.build(self.output()) as writers:
            for entry in self.data():
                if not entry.is_valid():
                    logger.info("Skipping invalid entry")
                self.write(writers, entry)


class Importer(EnsemblImporter):
    input_file = luigi.Parameter()
    destination = luigi.Parameter(default='/tmp')
    test = luigi.BoolParameter(default=False, significant=False)

    def is_pseudogene(self, feature):
        notes = feature.qualifiers.get('note', [])
        return any('pseudogene' in n for n in notes)

    def is_ncrna(self, feature):
        return feature.type in {'misc_RNA', 'ncRNA'}

    def assembly_info(self, feature):
        def as_exon(location):
            return {
                'primary_start': location.start + 1,
                'primary_end': int(location.end),
                'complement': location.strand == -1,
            }

        parts = [feature.location]
        if hasattr(feature.location, 'parts'):
            parts = feature.location.parts
        return [as_exon(l) for l in parts]

    def references(self, data):
        return [{
            'author': "Andrew Yates, Wasiu Akanni, M. Ridwan Amode, Daniel Barrell, Konstantinos Billis, Denise Carvalho-Silva, Carla Cummins, Peter Clapham, Stephen Fitzgerald, Laurent Gil1 Carlos Garcín Girón, Leo Gordon, Thibaut Hourlier, Sarah E. Hunt, Sophie H. Janacek, Nathan Johnson, Thomas Juettemann, Stephen Keenan, Ilias Lavidas, Fergal J. Martin, Thomas Maurel, William McLaren, Daniel N. Murphy, Rishi Nag, Michael Nuhn, Anne Parker, Mateus Patricio, Miguel Pignatelli, Matthew Rahtz, Harpreet Singh Riat, Daniel Sheppard, Kieron Taylor, Anja Thormann, Alessandro Vullo, Steven P. Wilder, Amonida Zadissa, Ewan Birney, Jennifer Harrow, Matthieu Muffato, Emily Perry, Magali Ruffier, Giulietta Spudich, Stephen J. Trevanion, Fiona Cunningham, Bronwen L. Aken, Daniel R. Zerbino, Paul Flicek",
            'location': "Nucleic Acids Res. 2016 44 Database issue:D710-6",
            'title': "Ensembl 2016",
            'pmid': 26687719,
            'doi': "10.1093/nar/gkv115",
        }]

    def standard_annotations(self, record):
        pattern = re.compile('\((.+)\)$')
        common_name = None
        species = record.annotations['organism']
        match = re.search(pattern, species)
        if match:
            common_name = match.group(1)
            species = re.sub(pattern, '', species).strip()

        taxid = None
        source = record.features[0]
        if source.type == 'source':
            taxid = int(qualifier_value(source, 'db_xref', '^taxon:(\d+)$'))

        return {
            'accession': None,
            'database': 'ensembl',
            'lineage': '; '.join(record.annotations['taxonomy']),
            'mol_type': 'genomic DNA',
            'pseudogene': 'N',
            'parent_accession': record.id,
            'seq_version': '',
            'common_name': common_name,
            'species': species,
            'ncbi_tax_id': taxid,
            'is_composite': 'N',
            'references': self.references(record),
            'sequence': record.seq,
        }

    def description(self, annotations, feature):
        species = annotations['species']
        if annotations.get('common_name', None):
            species += ' (%s)' % annotations['common_name']

        rna_type = self.ncrna(feature) or feature.type
        return '{species} {rna_type}'.format(
            species=species,
            rna_type=rna_type,
        )

    def gene(self, feature):
        return qualifier_value(feature, 'gene', '^(.+)$')

    def transcript(self, feature):
        return qualifier_value(feature, 'note', '^transcript_id=(.+)$')

    def ncrna(self, feature):
        return feature.qualifiers.get('note', [''])[0]

    def note(self, feature):
        note = feature.qualifiers.get('note', [])
        if len(note) > 1:
            return note[1:]
        return None

    def primary_id(self, annotations, feature):
        transcript = self.transcript(feature)
        ncrna = self.ncrna(feature)
        assert transcript, 'Cannot primary id without transcript'
        assert ncrna, 'Cannot primary id without ncRNA type'
        assert annotations['parent_accession']

        return '{parent}:{transcript}:{type}'.format(
            parent=annotations['parent_accession'],
            transcript=transcript,
            type=ncrna,
        )

    def sequence(self, annotations, feature):
        return str(feature.extract(annotations['sequence']))

    def entry_specific_data(self, feature):
        start, end = sorted([int(feature.location.start),
                             int(feature.location.end)])

        # Ensure that this is 1 based, even though biopython
        # switches to 0 based
        start += 1

        return {
            'primary_id': self.transcript(feature),
            'assembly_info': self.assembly_info(feature),
            'db_xrefs': feature.qualifiers.get('db_xref', []),
            'feature_location_start': start,
            'feature_location_end': end,
            'feature_type': feature.type,
            'gene': self.gene(feature),
            'ncrna_class': self.ncrna(feature),
            'note': self.note(feature),
            'locus_tag': feature.qualifiers.get('locus_tag', ''),
        }

    def rnacentral_entry(self, annotations, feature):
        data = dict(annotations)
        data.update(self.entry_specific_data(feature))
        data['accession'] = self.primary_id(annotations, feature)
        data['description'] = self.description(annotations, feature)
        data['sequence'] = self.sequence(annotations, feature)
        return RNAcentralEntry(**data)
