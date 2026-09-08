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

import json

import pytest
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.gtrnadb import helpers
from rnacentral_pipeline.databases.ncbi.taxonomy import TaxonomyEntry

# simple.json/version2.json are two single-entry extracts from the real
# GtRNAdb export (data/gtrnadb/other_eukaryotes_export_1.json), both for
# Aplysia californica (taxon 6500): a plain tRNA (simple) and one whose
# matureSequence differs from its genomic sequence (version2/"complex").


@pytest.fixture
def data():
    with open("data/gtrnadb/simple.json", "r") as raw:
        return json.load(raw)


@pytest.fixture
def data2():
    with open("data/gtrnadb/version2.json", "r") as raw:
        return json.load(raw)


@pytest.fixture
def taxonomy(tmp_path):
    db = SqliteDict(filename=str(tmp_path / "taxonomy.db"))
    db["6500"] = TaxonomyEntry(
        tax_id=6500,
        name="Aplysia californica",
        lineage=(
            "Eukaryota; Metazoa; Spiralia; Lophotrochozoa; Mollusca; "
            "Gastropoda; Heterobranchia; Euthyneura; Tectipleura; "
            "Aplysiida; Aplysioidea; Aplysiidae; Aplysia"
        ),
        aliases=[],
        replaced_by=None,
    )
    db.commit()
    return db


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


def test_url(data):
    assert (
        helpers.url(data[0])
        == "http://gtrnadb.ucsc.edu/genomes/eukaryota/Acali3/genes/tRNA-Ala-AGC-1-1.html"
    )


def test_anticodon(data):
    assert helpers.anticodon(data[0]) == "AGC"


def test_note_data(data):
    assert helpers.note_data(data[0]) == {
        "url": "http://gtrnadb.ucsc.edu/genomes/eukaryota/Acali3/genes/tRNA-Ala-AGC-1-1.html",
    }


def test_complex_note_data(data2):
    assert helpers.note_data(data2[0]) == {
        "url": "http://gtrnadb.ucsc.edu/genomes/eukaryota/Acali3/genes/tRNA-Arg-CCT-1-1.html",
    }


def test_lineage(data, taxonomy):
    assert helpers.lineage(taxonomy, data[0]) == (
        "Eukaryota; Metazoa; Spiralia; Lophotrochozoa; Mollusca; "
        "Gastropoda; Heterobranchia; Euthyneura; Tectipleura; "
        "Aplysiida; Aplysioidea; Aplysiidae; Aplysia"
    )


def test_species(data, taxonomy):
    assert helpers.species(taxonomy, data[0]) == "Aplysia californica"


def test_product(data):
    assert helpers.product(data[0]) == "tRNA-Ala (AGC)"


def test_as_dotbracket(data):
    ans = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))."
    assert helpers.dot_bracket(data[0]) == ans


def test_simple_description(data, taxonomy):
    assert (
        helpers.description(taxonomy, data[0]) == "Aplysia californica tRNA-Ala (AGC)"
    )


def test_complex_description(data2, taxonomy):
    assert (
        helpers.description(taxonomy, data2[0]) == "Aplysia californica tRNA-Arg (CCT)"
    )


def test_as_dotbracket_detects_weird_strings():
    data = {"secondaryStructure": ">>>...A<<<"}
    with pytest.raises(helpers.InvalidDotBracket):
        helpers.dot_bracket(data)


def test_primary_id_is_always_unique(data, data2):
    seen = set()
    possible = data + data2
    for entry in possible:
        for location in entry["genomeLocations"]:
            pid = helpers.primary_id(entry, location)
            assert pid not in seen
            seen.add(pid)
    assert seen


def test_builds_primary_id(data):
    pids = []
    entry = data[0]
    for location in entry["genomeLocations"]:
        pid = helpers.primary_id(entry, location)
        pids.append(pid)
    assert pids == ["GTRNADB:tRNA-Ala-AGC-1-1:KB942240.1:40028-40100"]


def test_chromosome(data):
    assert helpers.chromosome(data[0]["genomeLocations"][0]) == "scaffold00844"


def test_sequence_without_mature(data):
    assert (
        helpers.sequence(data[0])
        == "GGGGCTGTAGCTCAGGTGGTAGAGCGCTCGCTTAGCATGTGAGAGGTACCGGGATCGATACCCGGCAGCTCCA"
    )


def test_sequence_with_mature(data2):
    assert (
        helpers.sequence(data2[0])
        == "GCCTCCGTGGCCTAATGGATAAGGCATCGGCCTCCTAAGCCGGGGATTGCGGGTTCGAGTCCCGTCGGAGGTG"
    )
