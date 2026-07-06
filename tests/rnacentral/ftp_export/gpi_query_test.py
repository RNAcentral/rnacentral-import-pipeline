# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

from io import StringIO

from rnacentral_pipeline.rnacentral.ftp_export import gpi


def test_generic_query_no_filter():
    sql = str(gpi.generic_query(gpi.GpiFilter.NONE))
    assert "rnc_taxonomy" not in sql


def test_generic_query_species_filter():
    sql = str(gpi.generic_query(gpi.GpiFilter.SPECIES))
    assert '"rnc_taxonomy"' in sql
    assert "'species'" in sql


def test_generic_query_reference_proteome_filter():
    sql = str(gpi.generic_query(gpi.GpiFilter.REFERENCE_PROTEOME))
    assert '"rnc_taxonomy"' in sql
    assert "reference_proteome" in sql


def test_gpi_write_output_format():
    entries = [
        gpi.GpiEntry(
            urs_taxid="URS00001_9606",
            description="test description",
            rna_type="miRNA",
            symbol="test-sym",
            precursors=set(),
            aliases=[],
        ),
    ]
    out = StringIO()
    gpi.write(entries, out)
    lines = out.getvalue().splitlines()
    assert lines[0] == "!gpi-version: 1.2"
    fields = lines[1].split("\t")
    assert len(fields) == 10
    assert fields[0] == "RNAcentral"
    assert fields[1] == "URS00001_9606"
    assert fields[2] == "test-sym"
    assert fields[3] == "test description"
    assert fields[5] == "miRNA"
    assert fields[6] == "taxon:9606"
