import os
import re
import itertools as it

# from data import RNAcentralEntry

import attr
import luigi
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

from attr.validators import instance_of as is_a


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
class OutputWriters(object):
    short_sequences = attr.ib(validator=is_a(file))
    long_sequences = attr.ib(validator=is_a(file))
    references = attr.ib(validator=is_a(file))
    accessions = attr.ib(validator=is_a(file))
    locations = attr.ib(validator=is_a(file))

    @classmethod
    def build(cls, output):
        return cls(
            short_sequences=open(output.short_sequences, 'wb'),
            long_sequences=open(output.long_sequences, 'wb'),
            references=open(output.references, 'wb'),
            accessions=open(output.accesion, 'wb'),
            locations=open(output.locations, 'wb'),
        )

    def __enter__(self):
        return self

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
        with OutputWriters.build(self.output()) as writers:
            self.write(writers, self.data())


class Importer(EnsemblImporter):
    input_file = luigi.Parameter()
    destination = luigi.Parameter(default='/tmp')
    test = luigi.BoolParameter(default=False, significant=False)

    def is_pseudogene(self, feature):
        notes = feature.qualifiers.get('note', [])
        return any('pseudogene' in n for n in notes)

    def is_ncrna(self, feature):
        return feature.type in {'misc_RNA', 'ncRNA'}

    def assemblies(self, data):
        return []

    def references(self, data):
        return []

    def standard_annotations(self, record):
        pattern = re.compile('\((.+)\)$')
        common_name = None
        species = record.annotations['organism']
        match = re.search(pattern, species)
        if match:
            common_name = match.group(1)
            species = re.sub(pattern, '', species).strip()

        return {
            'database': 'ensembl',
            'lineage': '; '.join(record.annotations['taxonomy']),
            'mol_type': 'genomic DNA',
            'pseudogene': 'N',
            'parent_accession': record.id,
            'seq_version': '',
            'common_name': common_name,
            'species': species,
        }

    def description(self, anntotations, feature):
        return '$species $rna_type $product'.format(
            species=anntotations['species'],
            rna_type=None,
            product=None
        )

    def gene(self, feature):
        gene = feature.qualifiers.get('gene', [])
        if gene:
            return gene[0]
        return None

    def transcript(self, feature):
        transcripts = set()
        for note in feature.qualifiers.get('note', []):
            match = re.match('^transcript_id=(.+)$', note)
            if match:
                transcripts.add(match.group(1))
        if not transcripts:
            return None
        if len(transcripts) > 1:
            raise ValueError("Multiple transcripts for %s" % str(feature))
        return transcripts.pop()

    def ncrna(self, feature):
        return feature.qualifiers.get('note', [''])[0]

    def note(self, feature):
        note = feature.qualifiers.get('note', [])
        if len(note) > 1:
            return note[1:]
        return None

    def is_composite(self, feature):
        if issubclass(feature.location.__class__, CompoundLocation):
            return 'Y'
        return 'N'

    def create_rnacentral_entries(self, annotations, feature):
        return RNAcentralEntry.with_constants(
            annotations,
            primary_id=feature.id,
            sequence=str(feature.seq),
            ncbi_tax_id=annotations['taxid'],
            accession=None,
            assembly_info=self.assemblies(feature),
            db_xrefs=feature.qualifiers.get('db_xref', []),
            description=self.description(annotations, feature),
            feature_location_end=int(feature.location.end),
            feature_location_start=int(feature.location.start),
            feature_type=feature.type,
            gene=self.gene(feature),
            is_composite=self.is_composite(feature),
            ncrna_class=self.ncrna(feature),
            note=self.note(feature),
            optional_id='',
            references=self.references(feature),
            standard_name='',
            locus_tag=feature.qualifiers.get('locus_tag', ''),
        )
