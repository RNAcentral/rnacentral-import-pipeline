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

from databases import data
from databases.helpers import phylogeny as phy
from databases.helpers import publications as pub


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


def as_exon(assembly, exon):
    """
    Turn a raw exon into one we can store in the rnc_coordinates table.
    """

    complement = None
    if exon['strand'] == '+':
        complement = False
    elif exon['strand'] == '-':
        complement = True
    else:
        raise ValueError("Invalid strand %s" % exon)

    chromosome = exon['chromosome']
    if chromosome.startswith('chr'):
        chromosome = chromosome[3:]

    return data.Exon(
        chromosome_name=chromosome,
        # Input is 0 based, but we store 1 based
        primary_start=int(exon['startPosition']) + 1,
        primary_end=int(exon['endPosition']),
        assembly_id=assembly,
        complement=complement,
    )


def exons(genome_location):
    """
    Get all genomic locations this record is in.
    """
    assembly = genome_location['assembly']
    return [as_exon(assembly, e) for e in genome_location['exons']]


def gene_info(_):
    """
    Extract the gene level information from the given record if present.
    """
    pass


def as_reference(ref):
    """
    Turn a raw reference (just a pmid) into a reference we can import.
    """
    if ref.startswith('PMID:'):
        pmid = int(ref[5:])
        return pub.reference(pmid=pmid)
    return None


def references(record):
    """
    Get a list of all References in this record.
    """
    refs = it.imap(as_reference, record.get('publications', []))
    return list(it.ifilter(None, refs))


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
            return add_organism_preifx(ncrna, raw_gene['symbol'])

    raise ValueError("Could not create a name for %s" % ncrna)


def gene_synonyms(ncrna):
    """
    Find all gene synonyms, if they exist.
    """
    return ncrna.get('gene', {}).get('synonyms', [])


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


def parent_accession(location):
    """
    Find the parent accession. This will be the INSDC_accession of the first
    exon in this location.
    """

    if not location['exons']:
        return None

    first_exon = location['exons'][0]
    if 'INSDC_accession' not in first_exon:
        return None
    return first_exon['INSDC_accession'].split('.')[0]


def as_entry(database, p_accession, parsed_exons, record):
    """
    Generate an Entry to import based off the database, exons and raw record.
    """

    return data.Entry(
        primary_id=external_id(record),
        accession=record['primaryId'],
        ncbi_tax_id=taxid(record),
        database=database,
        sequence=record['sequence'],
        exons=parsed_exons,
        rna_type=record['soTermId'],
        url=record['url'],
        description=description(record),
        seq_version=record.get('version', '1'),
        parent_accession=p_accession,
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
        locations = record.get('genomeLocations', [])
        if not locations:
            yield as_entry(database, None, [], record)

        for location in locations:
            parsed_exons = exons(location)
            p_accession = parent_accession(location)
            yield as_entry(database, p_accession, parsed_exons, record)
