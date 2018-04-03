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

import collections as coll

from databases import data
from databases.helpers import phylogeny as phy


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


def as_location(location):
    exons = []
    for exon in location['exons']:
        complement = None
        if exon['strand'] == '+':
            complement = False
        elif exon['strand'] == '-':
            complement = True
        else:
            raise ValueError("Invalid strand %s" % exon)

        exons.append(data.Exon(
            chromosome=exon['chromosome'],
            primary_start=int(exon['startPosition']),
            primary_end=int(exon['endPosition']),
            complement=complement,
        ))
    return exons


def locations(record):
    """
    Get all genomic locations this record is in.
    """
    return [as_location(location) for location in record['genomeLocations']]


def gene_info(_):
    """
    Extract the gene level information from the given record if present.
    """
    pass


def as_reference(ref):
    """
    Turn a raw reference (just a pmid) into a reference we can import.
    """
    return data.Reference(
        pmid=ref['pubMedId'],
    )


def references(record):
    """
    Get a list of all References in this record.
    """
    return [as_reference(ref) for ref in record.get('publications', [])]


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
    return phy.species(taxid(ncrna))


def common_name(ncrna):
    return phy.common_name(taxid(ncrna))


def lineage(ncrna):
    return phy.lineage(taxid(ncrna))


def add_organism_preifx(ncrna, suffix):
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
        gene = ncrna['gene']
        if 'name' in gene:
            return add_organism_preifx(ncrna, gene['name'])
        if 'symbol' in gene:
            return add_organism_preifx(ncrna, gene['symbol'])

    raise ValueError("Could not create a name for %s" % ncrna)


def gene_synonyms(ncrna):
    return ncrna.get('gene', {}).get('synonyms', [])


def gene(ncrna):
    gene_id = ncrna.get('gene', {}).get('geneId', None)
    if gene_id:
        return gene_id.split(':', 1)[1]
    return None


def locus_tag(ncrna):
    return ncrna.get('gene', {}).get('locusTag', None)


def parent_accession(ncrna):
    locations = ncrna.get('genomeLocations', [])
    if not locations:
        return None

    first_exon = locations[0]['exons'][0]
    if 'INSDC_accession' not in first_exon:
        return None
    return first_exon['INSDC_accession'].split('.')[0]


def as_entry(database, exons, record):
    return data.Entry(
        primary_id=external_id(record),
        accession=record['primaryId'],
        ncbi_tax_id=taxid(record),
        database=database,
        sequence=record['sequence'],
        exons=exons,
        rna_type=record['soTermId'],
        url=record['url'],
        description=description(record),
        seq_version=record.get('version', '1'),
        parent_accession=parent_accession(record),
        xref_data=xrefs(record),
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
    Parses the given dict. This assumes the data is formatted according to
    version 1.0 (or equivalent) of the RNAcentral JSON schema.
    """

    database = raw['metaData']['dataProvider']
    for record in raw['data']:
        for exons in locations(record):
            yield as_entry(database, exons, record)
