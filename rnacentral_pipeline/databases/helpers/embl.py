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

import re
import collections as coll

from Bio import SeqIO

from rnacentral_pipeline.databases.data import Exon, Reference
import rnacentral_pipeline.databases.helpers.phylogeny as phy


class MissingTaxId(Exception):
    """
    This is  raised when an operation which should have a NCBI taxon id should,
    but does not.
    """
    pass


class MissingSource(Exception):
    """
    This is raised if we cannot extract the source feature from an EMBL record.
    """
    pass


def grouped_annotations(raw, split):
    """
    Parse a raw string into a dict. This will produce a key value mappign where
    the key is everything before the first split and values is everything
    after. The mapping values will be lists. This will correct RNACentral to
    RNAcentral.
    """

    parsed = coll.defaultdict(set)
    for entry in raw:
        if split not in entry:
            continue
        key, value = entry.split(split, 1)
        if key == 'RNACentral':
            key = 'RNAcentral'
        parsed[key].add(value)
    return {k: sorted(v) for k, v in parsed.items()}


def qualifier_value(feature, name, pattern, max_allowed=1):
    """
    This will parse the qualifier field defined by the given name for the given
    feature. This will extract all values matching the given regex pattern. If
    max allowed is 1 then only one distinct value is allowed and a single value
    will be returned. Otherwise all values in a set will be returned.
    """

    values = set()
    for note in feature.qualifiers.get(name, []):
        match = re.match(pattern, note)
        if match:
            values.add(match.group(1))
    if max_allowed is not None and len(values) > max_allowed:
        raise ValueError("Multiple values (%s) for %s in %s" %
                         (', '.join(sorted(values)), name, str(feature)))

    if len(values) == 0:
        return None
    if max_allowed == 1:
        return values.pop()
    return values


def qualifier_string(feature, name, separator=' '):
    """
    Extract the qualifier values and then join all results with the given
    separator, default ' '.
    """

    values = feature.qualifiers.get(name, [])
    if not values:
        return None
    return separator.join(values)


def source_feature(record):
    """
    Get the source feature for the record. This feature must be the first
    feature in the file.
    """

    source = record.features[0]
    if source.type != 'source':
        raise MissingSource("No source for: %s" % record)
    return source


def taxid(record):
    """
    Get the taxon id of the given record. This will pull the first feature,
    which must be of the 'source' type to do so.
    """

    source = source_feature(record)
    value = qualifier_value(source, 'db_xref', r'^taxon:(\d+)$')
    if value is None:
        return 32644  # Unclassified sequence
        # raise MissingTaxId("No taxid found for %s" % record)
    value = int(value)
    if value < 0:
        return 32644  # Unclassified sequence
    return value


def common_name(record):
    """
    Look up the common name of the record. This will use the taxid to pull an
    annotated common_name.
    """
    return phy.common_name(taxid(record))


def species(record):
    """
    Get the species from the record. This will use the taxid to pull the
    species, and will not use the annotated species.
    """
    return phy.species(taxid(record))


def organism(record):
    """
    Get the annotated organism in the record.
    """
    return record.qualifiers['organism']


def lineage(record):
    """
    Extract the taxon id and then query a remote server for the lineage for the
    given taxon id. Results are cached because it is common to constantly query
    with only the same few taxon ids.
    """
    return phy.lineage(taxid(record))


def project(record):
    """
    Get the annotated project xref if one exists.
    """

    for xref in record.dbxrefs:
        if 'Project' in xref:
            return xref.split(':', 1)[1]
    return None


def description(record):
    """
    Get the description of this record.
    """
    return record.description


def as_exon(record, location):
    """
    Turn a Biopython location object in a record to an Exon.
    """

    accession = record.annotations['accessions'][0]
    parts = accession.split(':')
    assembly_id = parts[1]
    chromosome_name = parts[2]
    return Exon(
        chromosome_name=chromosome_name,
        primary_start=location.start + 1,
        primary_end=int(location.end),
        assembly_id=assembly_id,
        complement=location.strand == -1,
    )


def exons(record, feature):
    """
    Turn all parts of the location into a list of exons to import.
    """

    parts = [feature.location]
    if hasattr(feature.location, 'parts'):
        parts = feature.location.parts
    return [as_exon(record, l) for l in parts]


def experiment(feature):
    """
    Lookup the annotated experiment information.
    """

    return qualifier_string(feature, 'experiment')


def inference(feature):
    """
    Look up the annotated inference information. THis will return a single ' '
    separate string of all inferences.
    """

    values = qualifier_value(feature, 'inference', r'^(.+)$', max_allowed=None)
    if values:
        return ' '.join(values)
    return None


def seq_version(record):
    """
    Get the sequence version given the record. If no sequence_version is
    annotated in the record then None is used.
    """

    return str(record.annotations.get('sequence_version', None))


def division(record):
    """
    get the division given a record.
    """

    return phy.division(taxid(record))


def standard_name(feature):
    """
    Get the standard name of feature.
    """
    return qualifier_value(feature, 'standard_name', '^(.+)$')


def gene(feature):
    """
    Get the gene this feature is a part of.
    """
    return qualifier_value(feature, 'gene', '^(.+)$')


def xref_data(feature):
    """
    Get a dict of the xref data for this feature. This will parse the db_xref
    qualifier to produce a key value mapping.
    """
    raw = feature.qualifiers.get('db_xref', [])
    return grouped_annotations(raw, ':')


def locus_tag(feature):
    """
    Get the locus tag of this feature. If none is present then the empty string
    is returned.
    """
    return feature.qualifiers.get('locus_tag', [None])[0]


def old_locus_tag(feature):
    """
    Get the old_locus_tag of this feature.
    """
    return feature.qualifiers.get('old_locus_tag', [None])[0]


def as_reference(ref):
    """
    Convert a biopython reference to one of our references. This will do some
    minor cleanups on the parsed reference to normalize some datatypes and
    remove bad title descriptions.
    """

    pmid = None
    if ref.pubmed_id:
        pmid = int(ref.pubmed_id)

    title = ref.title
    if title == ';':
        title = None

    return Reference(
        authors=ref.authors,
        location=ref.journal,
        title=title,
        pmid=pmid,
        doi=None,
    )


def references(record):
    """
    Turn all references for a record into a our type of references. This will
    assign all references to a single acecssion.
    """

    refs = record.annotations.get('references', [])
    return [as_reference(ref) for ref in refs]


def is_gene(feature):
    """
    Check if this feature is a gene
    """
    return feature.type == 'gene'


def transcripts(handle):
    """
    Fetch all the transcripts in the given EMBL file.
    """

    for record in SeqIO.parse(handle, 'embl'):
        current_gene = None
        for feature in record.features:
            if is_gene(feature):
                current_gene = feature
                continue

            assert gene(feature) != gene(current_gene)
            yield (record, current_gene, feature)


def sequence(record, feature):
    """
    Extract the sequence from the given record and feature.
    """
    return str(feature.extract(record.seq))
