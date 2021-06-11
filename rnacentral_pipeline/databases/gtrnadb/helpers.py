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
import urllib

import six

import requests

from rnacentral_pipeline.databases.data import Reference
import rnacentral_pipeline.databases.helpers.phylogeny as phy


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
    """
    Adds a http:// to the start of the url in data.
    """
    return "http://%s" % data["metadata"]["url"]


def anticodon(data):
    """
    Get the anticodon of this entry.
    """
    return data["metadata"]["anticodon"]


def note_data(data):
    """
    Create the dict that will be stored as a note. This is basically whatever
    GtRNAdb gives us, with a few duplicate or unneeded fields removed.
    """

    note = {}
    note.update(data["metadata"])
    extra = ["mature_sequence", "description"]
    for entry in extra:
        if entry in note:
            del note[entry]
    del note["organism"]
    del note["pseudogene"]
    for position in note["anticodon_positions"]:
        position["relative_start"] = int(position["relative_start"])
        position["relative_stop"] = int(position["relative_stop"])
    note["score"] = float(note["score"])
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


def lineage(taxonomy, data):
    """
    Get a standardized lineage for the given taxon id.
    """
    if data["ncbi_tax_id"] in taxonomy:
        return taxonomy["ncbi_tax_id"].lineage
    return phy.lineage(data["ncbi_tax_id"])


def species(taxonomy, data):
    """
    Get a standardized species name for the given taxon id.
    """
    if data["ncbi_tax_id"] in taxonomy:
        return taxonomy["ncbi_tax_id"].name
    return phy.species(data["ncbi_tax_id"])


def description(taxonomy, data):
    """
    Generate a description for the entries specified by the data.
    """
    details = data["metadata"].get("description", product(data))
    return "{name} {details}".format(
        name=species(taxonomy, data),
        details=details,
    )


def product(data):
    """
    Generate the product for the entries specified by the data.
    """
    return "tRNA-{aa} ({anticodon})".format(
        aa=data["metadata"]["isotype"],
        anticodon=data["metadata"]["anticodon"],
    )


def primary_id(data, location):
    """
    Generate a primary key for the given data and location.
    """

    start = min(int(e["start"]) for e in location["exons"])
    stop = max(int(e["stop"]) for e in location["exons"])
    return "{gene}:{accession}:{start}-{stop}".format(
        gene=data["gene"],
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
        gene=data["gene"],
    )


def seq_version(_):
    """
    Compute a seq_version for GtRNAdb data. CUrrentlyt his will always return
    '1'
    """
    return "1"


def references():
    """
    Returns the default accessions for GtRNAdb data.
    """

    return [
        Reference(
            authors="Chan P.P., Lowe T.M.",
            location="Nucl. Acids Res. 37(Database issue)",
            title=(
                "GtRNAdb: A database of transfer RNA genes detected in "
                "genomic sequence"
            ),
            pmid=18984615,
            doi=u"10.1093/nar/gkn787.",
        )
    ]


def sequence(data):
    if "mature_sequence" in data["metadata"]:
        return data["metadata"]["mature_sequence"].upper()
    return data["sequence"].upper()


def features(data):
    return []


def regions(data):
    pass
