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

import collections as coll
import re
import typing as ty

from Bio import SeqIO

import rnacentral_pipeline.databases.helpers.phylogeny as phy
import rnacentral_pipeline.databases.helpers.publications as pubs
from rnacentral_pipeline.databases import data

IGNORE_FEATURES = {
    "source",
    "exon",
    "STS",
    "misc_feature",
}


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


def grouped_annotations(raw, split) -> ty.Dict[str, ty.List[str]]:
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
        if key == "RNACentral":
            key = "RNAcentral"
        parsed[key].add(value)
    return {k: sorted(v) for k, v in parsed.items()}


def qualifier_value(feature, name, pattern, max_allowed=1) -> ty.Optional[str]:
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
        raise ValueError(
            "Multiple values (%s) for %s in %s"
            % (", ".join(sorted(values)), name, str(feature))
        )

    if len(values) == 0:
        return None
    if max_allowed == 1:
        return values.pop()
    return values


def qualifier_string(feature, name, separator=" ") -> ty.Optional[str]:
    """
    Extract the qualifier values and then join all results with the given
    separator, default ' '.
    """

    values = feature.qualifiers.get(name, [])
    if not values:
        return None
    return separator.join(values)


def source_qualifier_value(
    record, qualifier, pattern=r"^(.+)$", **kwargs
) -> ty.List[str]:
    source = source_feature(record)
    return qualifier_value(source, qualifier, pattern, **kwargs)


def source_feature(record):
    """
    Get the source feature for the record. This feature must be the first
    feature in the file.
    """

    source = record.features[0]
    if source.type != "source":
        raise MissingSource("No source for: %s" % record)
    return source


def taxid(record) -> int:
    """
    Get the taxon id of the given record. This will pull the first feature,
    which must be of the 'source' type to do so.
    """

    source = source_feature(record)
    value = qualifier_value(source, "db_xref", r"^taxon:(\d+)$")
    if value is None:
        return 32644  # Unclassified sequence
        # raise MissingTaxId("No taxid found for %s" % record)
    value = int(value)
    if value < 0:
        return 32644  # Unclassified sequence
    return value


def common_name(record) -> str:
    """
    Look up the common name of the record. This will use the taxid to pull an
    annotated common_name.
    """
    return phy.common_name(taxid(record))


def species(record) -> str:
    """
    Get the species from the record. This will use the taxid to pull the
    species, and will not use the annotated species.
    """
    return phy.species(taxid(record))


def organism(record) -> str:
    """
    Get the annotated organism in the record.
    """
    return record.qualifiers["organism"]


def lineage(record) -> str:
    """
    Extract the taxon id and then query a remote server for the lineage for the
    given taxon id. Results are cached because it is common to constantly query
    with only the same few taxon ids.
    """
    return phy.lineage(taxid(record))


def project(record) -> ty.Optional[str]:
    """
    Get the annotated project xref if one exists.
    """

    for xref in record.dbxrefs:
        if "Project" in xref:
            return xref.split(":", 1)[1]
    return None


def description(record) -> str:
    """
    Get the description of this record.
    """
    return record.description


def experiment(feature) -> ty.Optional[str]:
    """
    Lookup the annotated experiment information.
    """

    return qualifier_string(feature, "experiment")


def inference(feature) -> ty.Optional[str]:
    """
    Look up the annotated inference information. THis will return a single ' '
    separate string of all inferences.
    """

    values = qualifier_value(feature, "inference", r"^(.+)$", max_allowed=None)
    if values:
        return " ".join(values)
    return None


def seq_version(record):
    """
    Get the sequence version given the record. If no sequence_version is
    annotated in the record then None is used.
    """

    return str(record.annotations.get("sequence_version", None))


def division(record):
    """
    get the division given a record.
    """

    return phy.division(taxid(record))


def standard_name(feature):
    """
    Get the standard name of feature.
    """
    return qualifier_value(feature, "standard_name", "^(.+)$")


def gene(feature):
    """
    Get the gene this feature is a part of.
    """
    return qualifier_value(feature, "gene", "^(.+)$")


def xref_data(feature):
    """
    Get a dict of the xref data for this feature. This will parse the db_xref
    qualifier to produce a key value mapping.
    """
    raw = feature.qualifiers.get("db_xref", [])
    return grouped_annotations(raw, ":")


def locus_tag(feature):
    """
    Get the locus tag of this feature. If none is present then the empty string
    is returned.
    """
    return feature.qualifiers.get("locus_tag", [None])[0]


def as_reference(ref):
    """
    Convert a biopython reference to one of our references. This will do some
    minor cleanups on the parsed reference to normalize some datatypes and
    remove bad title descriptions.
    """

    pmid = None
    if ref.pubmed_id:
        return pubs.reference(int(ref.pubmed_id))

    title = ref.title
    if title == ";":
        title = None
    else:
        title = str(title)

    return data.Reference(
        authors=str(ref.authors),
        location=str(ref.journal),
        title=title,
        pmid=pmid,
        doi=None,
    )


def references(record):
    """
    Turn all references for a record into a our type of references. This will
    assign all references to a single acecssion.
    """

    refs = record.annotations.get("references", [])
    return [as_reference(ref) for ref in refs]


def is_gene(feature):
    """
    Check if this feature is a gene
    """
    return feature.type == "gene"


def rna_type(feature):
    """
    Compute the RNA type of the given feature. This ensures that the feature
    has a known ncRNA feature type.
    """

    if feature.type == "ncRNA":
        return qualifier_value(feature, "ncRNA_class", r"^(.+)$")
    elif feature.type in {"misc_RNA", "rRNA", "tRNA", "precursor_RNA"}:
        return feature.type

    raise ValueError("Non-ncRNA feature type: %s" % feature.type)


def transcripts(handle):
    """
    Fetch all the transcripts in the given EMBL file.
    """

    for record in SeqIO.parse(handle, "embl"):
        current_gene = None
        for feature in record.features:
            if feature.type in IGNORE_FEATURES:
                continue

            if is_gene(feature):
                current_gene = feature
                continue

            assert current_gene is None or gene(feature) == gene(current_gene)
            yield (record, current_gene, feature)


def sequence(record, feature):
    """
    Extract the sequence from the given record and feature.
    """
    return str(feature.extract(record.seq))


def mol_type(source):
    return qualifier_value(source, "mol_type", r"^(.+)$")


def chromosome(source):
    chromosomes = source.qualifiers["chromosome"]
    return chromosomes[0]


def gene_synonyms(feature):
    result = []
    synonyms = feature.qualifiers.get("gene_synonym", [])
    for synonym in synonyms:
        result.extend(synonym.split("; "))
    return result


def product(feature):
    return qualifier_string(feature, "product", separator="; ")


def organelle(source):
    values = qualifier_value(source, "organelle", r"^(.+)$", max_allowed=None)
    if not values:
        return None
    if len(values) == 1:
        return values.pop()
    if len(values) == 2:
        # Seems strange but there are cases where the value is
        # ['plastid', 'plastid:chloroplast'] and in that case we want
        # 'plastid:chloroplast' as that is more specific.
        first, second = sorted(values, key=len)
        if second.startswith(first):
            return second
    return " ".join(sorted(values))
