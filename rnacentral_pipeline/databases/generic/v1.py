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

import itertools as it
import collections as coll

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pub


class UnexpectedCoordinates(Exception):
    """
    Raised if we get coordinates to store but we do not know what the coordinate
    system used is.
    """
    pass


@attr.s()
class Context(object):
    database = attr.ib(validator=is_a(str))
    coordinate_system = attr.ib(validator=optional(is_a(data.CoordinateSystem)))


def secondary_structure(record):
    """
    Fetches the secondary structure, if any, of the given JSON schema entry.
    """

    dot_bracket = record.get('secondary_structure', None)
    if dot_bracket:
        return data.SecondaryStructure(dot_bracket=dot_bracket)
    return data.SecondaryStructure.empty()


def xrefs(record):
    """
    Fetches the cross references between this sequence and any other known
    databases.
    """

    grouped = coll.defaultdict(list)
    for xref in record.get('crossReferenceIds', []):
        dbid, accession = xref.split(':', 1)
        grouped[dbid].append(accession)
    return dict(grouped)


def taxid(entry):
    """
    Gets the NCBI taxon id as an integer.
    """

    base, tid = entry['taxonId'].split(':', 1)
    assert base == 'NCBITaxon'
    return int(tid)


def as_exon(exon, context):
    """
    Turn a raw exon into one we can store in the rnc_coordinates table.
    """

    start_fix = 0
    stop_fix = 0
    if context.database == 'MIRBASE':
        start_fix = 1
        stop_fix = 1

    return data.Exon(
        start=int(exon['startPosition']) + start_fix,
        stop=int(exon['endPosition']) + stop_fix,
    )


def as_region(region, context):
    """
    Turn a raw region in the JSON document into a SequenceRegion object.
    """

    exons = region['exons']
    chromosome = exons[0]['chromosome']
    if chromosome.startswith('chr'):
        chromosome = chromosome[3:]

    return data.SequenceRegion(
        chromosome=chromosome,
        strand=exons[0]['strand'],
        exons=[as_exon(e, context) for e in exons],
        assembly_id=region['assembly'],
        coordinate_system=context.coordinate_system,
    )


def regions(entry, context):
    """
    Get all genomic locations this record is in.
    """

    result = []
    for region in entry.get('genomeLocations', []):
        if not region.get('exons', None):
            continue

        if not context.coordinate_system:
            raise UnexpectedCoordinates(region)

        result.append(as_region(region, context))
    return result


def gene_info(_):
    """
    Extract the gene level information from the given record if present.
    """
    pass


def references(record):
    """
    Get a list of all References in this record.
    """
    return [pub.reference(r) for r in record.get('publications', [])]


def anticodon(record):
    """
    Get the anticodon information, if any, of this record.
    """
    return record.get('sequenceFeatures', {}).get('anticodon', None)


def external_id(record):
    """
    Extract the external id for the record.
    """

    _, eid = record['primaryId'].split(':', 1)
    return eid


def species(ncrna):
    """
    Get the specis information using the taxid.
    """
    return phy.species(taxid(ncrna))


def common_name(ncrna):
    """
    Get the common name using the taxid.
    """
    return phy.common_name(taxid(ncrna))


def lineage(ncrna):
    """
    Look up the lineage using the taxid.
    """
    return phy.lineage(taxid(ncrna))


def add_organism_preifx(ncrna, suffix):
    """
    Add a the 'species (common name)' prefix to a generic suffix. This is
    commonly used to build nice descriptions.
    """

    prefix = species(ncrna)
    name = common_name(ncrna)
    if name:
        prefix += ' (%s)' % name
    return prefix + ' ' + suffix


def description(ncrna):
    """
    Generate a description for the given ncrna record. This will use a series
    of choices to select a good description.
    """

    if 'description' in ncrna and ncrna['description']:
        return ncrna['description']

    if 'name' in ncrna and ncrna['name']:
        return add_organism_preifx(ncrna, ncrna['name'])

    if 'gene' in ncrna:
        raw_gene = ncrna['gene']
        if 'name' in raw_gene:
            return add_organism_preifx(ncrna, raw_gene['name'])
        if 'symbol' in raw_gene:
            symbol = raw_gene['symbol']
            db_prefix = ncrna['primaryId'].split(':', 1)[0]
            if symbol.startswith(db_prefix):
                symbol = symbol.split(':', 1)[1]
            return add_organism_preifx(ncrna, symbol)

    raise ValueError("Could not create a name for %s" % ncrna)


def gene_synonyms(ncrna):
    """
    Find all gene synonyms, if they exist.
    """
    gene = ncrna.get('gene', {})
    synonyms = gene.get('synonyms', [])
    if 'symbol' in gene:
        synonyms.append(gene['symbol'])
    return synonyms


def gene(ncrna):
    """
    Get the id of the gene, if possible.
    """

    gene_id = ncrna.get('gene', {}).get('geneId', None)
    if gene_id:
        return gene_id.split(':', 1)[1]
    return None


def locus_tag(ncrna):
    """
    Get the locus tag of the gene, if present.
    """
    return ncrna.get('gene', {}).get('locusTag', None)


def optional_id(record, context):
    """
    Create an optional id for mirbase entries. This basically uses the name
    field, which will be the miRBase gene name.
    """

    if 'description' in record and \
            'name' in record and ' ' not in record['name']:
        return record['name']
    if context.database == 'MIRBASE':
        return record['name']
    return None


def related_sequences(record):
    """
    This will create a list of RelatedSequences from the given record.
    """

    sequences = []
    for related in record.get('relatedSequences', []):
        evidence = related.get('evidence', {})
        if evidence:
            # pylint: disable=star-args
            evidence = data.RelatedEvidence(**evidence)
        else:
            evidence = data.RelatedEvidence.empty()

        coordinates = []
        for coord in related.get('coordinates', []):
            coordinates.append(data.RelatedCoordinate(
                start=coord['startPosition'],
                stop=coord['endPosition'],
            ))

        sequences.append(data.RelatedSequence(
            sequence_id=related['sequenceId'],
            relationship=related['relationship'],
            coordinates=coordinates,
            evidence=evidence,
        ))
    return sequences


def add_related_by_gene(entries):
    """
    This will modify the related sequences so that it contains between all
    transcripts that have the same gene id.
    """

    updated = []
    for first in entries:
        related = []
        for second in entries:
            if first.accession == second.accession:
                continue

            related.append(data.RelatedSequence(
                sequence_id=second.accession,
                relationship='isoform',
            ))

        updated.append(attr.evolve(
            first,
            related_sequences=first.related_sequences + related,
        ))
    return updated


def note_data(record):
    return {
        'url': record['url'],
    }


def coordinate_system(metadata):
    system = metadata.get('genomicCoordinateSystem', '1-start, fully-closed')
    if not system:
        raise ValueError("Could not find coordinate system for: %s" % metadata)
    return data.CoordinateSystem.from_name(system)


def as_entry(record, context):
    """
    Generate an Entry to import based off the database, exons and raw record.
    """
    return data.Entry(
        primary_id=external_id(record),
        accession=record['primaryId'],
        ncbi_tax_id=taxid(record),
        database=context.database,
        sequence=record['sequence'],
        regions=regions(record, context),
        rna_type=record['soTermId'],
        url=record['url'],
        seq_version=record.get('version', '1'),
        optional_id=optional_id(record, context),
        description=description(record),
        note_data=note_data(record),
        xref_data=xrefs(record),
        related_sequences=related_sequences(record),
        species=species(record),
        lineage=lineage(record),
        common_name=common_name(record),
        secondary_structure=secondary_structure(record),
        references=references(record),
        organelle=record.get('localization', None),
        product=record.get('product', None),
        anticodon=anticodon(record),
        gene=gene(record),
        gene_synonyms=gene_synonyms(record),
        locus_tag=locus_tag(record),
    )


def parse(raw):
    """
    Parses the given dict into a Entry object. This assumes the data is
    formatted according to version 1.0 (or equivalent) of the RNAcentral JSON
    schema.
    """

    def key(raw):
        return gene(raw) or ''

    context = Context(
        database=raw['metaData']['dataProvider'],
        coordinate_system=coordinate_system(raw['metaData']),
    )

    ncrnas = sorted(raw['data'], key=key)

    metadata_pubs = raw['metaData'].get('publications', [])
    metadata_refs = [pub.reference(r) for r in metadata_pubs]

    for gene_id, records in it.groupby(ncrnas, gene):
        entries = [as_entry(r, context) for r in records]

        if gene_id:
            entries = add_related_by_gene(entries)

        for entry in entries:
            refs = entry.references + metadata_refs
            yield attr.evolve(entry, references=refs)
