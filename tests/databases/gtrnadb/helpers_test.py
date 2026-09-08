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

"""
GtRNAdb changed schema in 2021 (c5b3d9c0) and the old simple.json/version2.json
fixtures went with it, which left most of this file loading files that no longer
exist. These are the same helpers rewritten against the export we actually
parse, with the expectations read off the record rather than off the code.
"""

import json

import pytest

import rnacentral_pipeline.databases.helpers.publications as pub
from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.gtrnadb import helpers

EXPORT = "data/gtrnadb/other_eukaryotes_export_1.json"


@pytest.fixture(scope="module")
def export():
    with open(EXPORT, "r") as raw:
        return json.load(raw)


@pytest.fixture(scope="module")
def records(export):
    return export["data"]


@pytest.fixture
def record(records):
    """A single exon tRNA whose mature and genomic sequences agree."""
    return next(r for r in records if r["primaryId"] == "GTRNADB:tRNA-Ala-AGC-1-1")


@pytest.fixture
def spliced(records):
    """A two exon tRNA, so the intron is spliced out of matureSequence."""
    return next(r for r in records if r["primaryId"] == "GTRNADB:tRNA-Arg-CCT-1-1")


def test_can_find_all_remote_urls():
    assert (
        helpers.extract_download_urls(
            "http://google.com",
            """
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
 <head>
  <title>Index of /download/RNAcentral</title>
 </head>
 <body>
<h1>Index of /download/RNAcentral</h1>
<pre><img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                    <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[DIR]"> <a href="/download/">Parent Directory</a>                             -
<img src="/icons/compressed.gif" alt="[   ]"> <a href="archaea_tRNAs.json.gz">archaea_tRNAs.json.gz</a>   21-Nov-2017 08:24  934K
<img src="/icons/compressed.gif" alt="[   ]"> <a href="bacteria_tRNAs.tar.gz">bacteria_tRNAs.tar.gz</a>   22-Aug-2017 02:29   15M
<img src="/icons/compressed.gif" alt="[   ]"> <a href="fungi_tRNAs.tar.gz">fungi_tRNAs.tar.gz</a>      21-Nov-2017 08:25  5.7M
<img src="/icons/compressed.gif" alt="[   ]"> <a href="model_tRNAs.tar.gz">model_tRNAs.tar.gz</a>      24-Nov-2017 00:31  126K
<hr></pre>
<address>Apache/2.2.15 (CentOS) Server at <a href="mailto:lowe@soe.ucsc.edu">trna.ucsc.edu</a> Port 80</address>
</body></html>
    """,
        )
        == [
            ("archaea_tRNAs.json.gz", "http://google.com/archaea_tRNAs.json.gz"),
            ("bacteria_tRNAs.tar.gz", "http://google.com/bacteria_tRNAs.tar.gz"),
            ("fungi_tRNAs.tar.gz", "http://google.com/fungi_tRNAs.tar.gz"),
            ("model_tRNAs.tar.gz", "http://google.com/model_tRNAs.tar.gz"),
        ]
    )


def test_complains_if_no_download_urls():
    with pytest.raises(Exception):
        assert helpers.extract_download_urls(
            "http://google.com",
            """
    <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
    <html>
     <head>
      <title>Index of /download/RNAcentral</title>
     </head>
     <body>
    <h1>Index of /download/RNAcentral</h1>
    <pre><img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>                    <a href="?C=M;O=A">Last modified</a>      <a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr><img src="/icons/back.gif" alt="[DIR]"> <a href="/download/">Parent Directory</a>                             -
    <img src="/icons/compressed.gif" alt="[   ]"> <a >archaea_tRNAs.json.gz</a>   21-Nov-2017 08:24  934K
    <img src="/icons/compressed.gif" alt="[   ]"> <a >bacteria_tRNAs.tar.gz</a>   22-Aug-2017 02:29   15M
    <img src="/icons/compressed.gif" alt="[   ]"> <a >fungi_tRNAs.tar.gz</a>      21-Nov-2017 08:25  5.7M
    <img src="/icons/compressed.gif" alt="[   ]"> <a >model_tRNAs.tar.gz</a>      24-Nov-2017 00:31  126K
    <hr></pre>
    <address>Apache/2.2.15 (CentOS) Server at <a href="mailto:lowe@soe.ucsc.edu">trna.ucsc.edu</a> Port 80</address>
    </body></html>
        """,
        )


def test_url(record):
    assert helpers.url(record) == (
        "http://gtrnadb.ucsc.edu/genomes/eukaryota/Acali3/genes/tRNA-Ala-AGC-1-1.html"
    )


def test_anticodon(record):
    assert helpers.anticodon(record) == "AGC"


def test_note_data_is_just_the_url(record):
    assert helpers.note_data(record) == {"url": helpers.url(record)}


def test_taxid(record):
    assert helpers.taxid(record) == 6500


def test_product(record):
    assert helpers.product(record) == "tRNA-Ala (AGC)"


def test_gene_naming(record):
    assert helpers.gene(record) == "tRNA-Ala-AGC-1-1"
    assert helpers.optional_id(record) == "tRNA-Ala-AGC-1-1"
    assert helpers.gene_synonyms(record) == ["scaffold00844.trna1-AlaAGC"]


def test_chromosome(record):
    assert helpers.chromosome(record["genomeLocations"][0]) == "scaffold00844"


def test_chromosome_renames_the_bare_name():
    location = {"exons": [{"chromosome": "Chromosome"}]}
    assert helpers.chromosome(location) == "chr"


def test_accessions(record):
    location = record["genomeLocations"][0]
    assert helpers.parent_accession(location) == "KB942240.1"
    assert helpers.accession(record, location) == "KB942240.1:tRNA-Ala-AGC-1-1"
    assert helpers.seq_version(record) == "1"


def test_builds_primary_id(record):
    location = record["genomeLocations"][0]
    assert helpers.primary_id(record, location) == (
        "GTRNADB:tRNA-Ala-AGC-1-1:KB942240.1:40028-40100"
    )


def test_primary_id_spans_every_exon(spliced):
    location = spliced["genomeLocations"][0]
    # Exons at 1827195-1827231 and 1827257-1827292, so the id covers both.
    assert helpers.primary_id(spliced, location) == (
        "GTRNADB:tRNA-Arg-CCT-1-1:KB941428.1:1827195-1827292"
    )


def test_primary_id_is_always_unique(records):
    seen = set()
    for entry in records:
        for location in entry["genomeLocations"]:
            pid = helpers.primary_id(entry, location)
            assert pid not in seen
            seen.add(pid)
    assert seen


def test_sequence_prefers_the_mature_form(spliced):
    assert spliced["matureSequence"] != spliced["sequence"]
    assert helpers.sequence(spliced) == spliced["matureSequence"].upper()


def test_sequence_falls_back_to_the_genomic_form(record):
    without_mature = {k: v for k, v in record.items() if k != "matureSequence"}
    assert helpers.sequence(without_mature) == record["sequence"].upper()


def test_features_carry_the_anticodon(record):
    assert helpers.features(record) == [
        data.SequenceFeature(
            name="anticodon",
            feature_type="anticodon",
            location=[34, 35, 36],
            sequence="AGC",
            provider="GTRNADB",
            metadata={"isotype": "Ala", "sequence": "AGC"},
        )
    ]


def test_no_features_without_an_anticodon(record):
    stripped = dict(record, sequenceFeatures={"isotype": "Ala"})
    assert helpers.features(stripped) == []


def test_regions_are_one_based(record):
    region = helpers.regions(record["genomeLocations"][0])[0]
    assert region.assembly_id == "AplCal3.0"
    assert region.chromosome == "scaffold00844"
    assert region.strand == data.Strand.forward
    assert [(e.start, e.stop) for e in region.exons] == [(40028, 40100)]


def test_regions_keeps_every_exon(spliced):
    region = helpers.regions(spliced["genomeLocations"][0])[0]
    assert [(e.start, e.stop) for e in region.exons] == [
        (1827195, 1827231),
        (1827257, 1827292),
    ]


def test_references_come_from_the_metadata(export):
    assert helpers.references(export["metaData"]) == [
        pub.reference("PMID:18984615"),
        pub.reference("PMID:26673694"),
    ]


def test_dot_bracket_translates_the_gtrnadb_notation():
    raw = {"secondary_structure": ">>>>>>>..>>>>........<<<<."}
    assert helpers.dot_bracket(raw) == "(((((((..((((........))))."


def test_dot_bracket_detects_weird_strings():
    with pytest.raises(helpers.InvalidDotBracket):
        helpers.dot_bracket({"secondary_structure": ">>>...A<<<"})
