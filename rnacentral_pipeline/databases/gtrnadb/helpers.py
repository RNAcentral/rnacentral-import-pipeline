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
import typing as ty

import six

import requests
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases import data
import rnacentral_pipeline.databases.helpers.phylogeny as phy
import rnacentral_pipeline.databases.helpers.publications as pub
from rnacentral_pipeline.databases.ncbi.taxonomy import TaxonomyEntry


class InvalidDotBracket(Exception):
    """
    This is raised when the given string cannot be turned into a valid
    dot-bracket string.
    """

    pass


def extract_download_urls(base_url, text):
    """
    Given a chunk of text returned by the main GtRNAdb RNAcentral download url
    (http://trna.ucsc.edu/download/RNAcentral/) this can pull out the URLs that
    represent the filenames to download.
    """

    data = []
    for line in text.split("\n"):
        # I am aware of how awful it is to parse HTML via regex, but I long for
        # the return of the elder gods and I'm doing my part to bring them
        # back.
        if "/icons/compressed.gif" in line:
            match = re.search('href="([^"]+)"', line)
            assert match
            filename = match.group(1)
            data.append((filename, six.moves.urllib_parse.urljoin(base_url, filename)))
    assert data, "Given text contained now downloadable urls"
    return data


def downloadable_files(base_url):
    """
    Determine all remote files to download.
    """

    response = requests.get(base_url)
    return extract_download_urls(base_url, response.text)


def url(data):
    return data["url"]


def anticodon(data):
    """
    Get the anticodon of this entry.
    """
    return data["sequenceFeatures"]["anticodon"]["sequence"]


def note_data(data):
    """
    Create the dict that will be stored as a note. This is basically whatever
    GtRNAdb gives us, with a few duplicate or unneeded fields removed.
    """

    note = {}
    note["url"] = url(data)
    return note


def chromosome(location):
    """
    Get the chromosome this location is part of.
    """

    chrom = location["exons"][0]["chromosome"]
    if chrom == "Chromosome":
        return "chr"
    return chrom


def taxid(data: ty.Dict[str, str]) -> int:
    _, taxid = data["taxonId"].split(":", 1)
    return int(taxid)


def lineage(taxonomy: SqliteDict, data):
    """
    Get a standardized lineage for the given taxon id.
    """
    ncbi_taxid = taxid(data)
    if str(ncbi_taxid) in taxonomy:
        return taxonomy[ncbi_taxid].lineage
    return phy.lineage(ncbi_taxid)


def species(taxonomy, data):
    """
    Get a standardized species name for the given taxon id.
    """
    ncbi_taxid = taxid(data)
    if ncbi_taxid in taxonomy:
        return taxonomy[ncbi_taxid].name
    return phy.species(ncbi_taxid)


def description(taxonomy, data):
    """
    Generate a description for the entries specified by the data.
    """
    return "{name} {details}".format(
        name=species(taxonomy, data),
        details=product(data),
    )


def product(data):
    """
    Generate the product for the entries specified by the data.
    """
    return data["name"]


def primary_id(data, location):
    """
    Generate a primary key for the given data and location.
    """

    start = min(int(e["startPosition"]) for e in location["exons"])
    stop = max(int(e["endPosition"]) for e in location["exons"])
    return "GTRNADB:{gene}:{accession}:{start}-{stop}".format(
        gene=data["gene"]["symbol"],
        accession=location["exons"][0]["INSDC_accession"],
        start=start,
        stop=stop,
    )


def dot_bracket(data):
    """
    Generate a dot bracket string from the secondary structure string that
    GtRNAdb uses. That is turn '>>..<<' to '((..))'.
    """

    transformed = data["secondary_structure"].replace(">", "(").replace("<", ")")

    if set(transformed) != set("(.)"):
        raise InvalidDotBracket("Unexpected characters in %s" % transformed)

    return transformed


def parent_accession(location):
    """
    Get the parent accessino for the given location.
    """
    return location["exons"][0]["INSDC_accession"]


def accession(data, location):
    """
    Generate an accession for the given location in data.
    """
    return "{ac}:{gene}".format(
        ac=parent_accession(location),
        gene=data["gene"]["symbol"],
    )


def seq_version(_):
    """
    Compute a seq_version for GtRNAdb data. CUrrentlyt his will always return
    '1'
    """
    return "1"


def references(metadata):
    """
    Returns the default accessions for GtRNAdb data.
    """

    return [pub.reference(pmid) for pmid in metadata["publications"]]


def sequence(data):
    if "matureSequence" in data:
        return data["matureSequence"].upper()
    return data["sequence"].upper()


def features(raw):
    anti = raw["sequenceFeatures"].get("anticodon", None)
    if not anti:
        return []

    return [
        data.SequenceFeature(
            name="anticodon",
            feature_type="anticodon",
            location=anti["indexes"],
            sequence=anti["sequence"],
            metadata={
                "isotype": raw["sequenceFeatures"]["isotype"],
                "sequence": anti["sequence"],
            },
        )
    ]


def regions(location):
    if not location:
        return []
    strand = data.Strand.build(location["exons"][0]["strand"])
    exons = []
    for exon in location["exons"]:
        start = exon["startPosition"]
        stop = exon["endPosition"]
        if start > stop:
            start, stop = stop, start
        exons.append(
            data.Exon(
                start=start,
                stop=stop,
            )
        )
    exons = tuple(exons)
    return [
        data.SequenceRegion(
            assembly_id=location["assembly"],
            chromosome=chromosome(location),
            strand=strand,
            exons=exons,
            coordinate_system=data.CoordinateSystem.one_based(),
        )
    ]


def gene_synonyms(raw: ty.Dict[str, ty.Any]) -> ty.List[str]:
    return raw["gene"]["synonyms"]


def gene(raw: ty.Dict[str, ty.Any]) -> str:
    return raw["gene"]["symbol"]


def optional_id(raw) -> str:
    return raw["gene"]["symbol"]
