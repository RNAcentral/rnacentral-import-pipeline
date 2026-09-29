# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

find-sequences.sql builds the set of URS whose types R2DT has no template for,
then selects what to draw. The upi -> urs rename (3e0b95e0) silently restored
an older final query that ignored that set, so lncRNAs went back to being
force-fitted and the per-row xref probe that once stalled extraction for 14h
came back.
"""

from pathlib import Path

SQL = (
    Path(__file__).resolve().parents[3] / "files/r2dt/find-sequences.sql"
).read_text()
EXCLUDED, FINAL = SQL.split("COPY (", 1)


def test_the_final_query_skips_the_excluded_set():
    assert "excluded" in FINAL


def test_the_final_query_does_not_probe_xref_per_sequence():
    assert "xref" not in FINAL


def test_whole_mrnas_are_never_sent_to_r2dt():
    assert "'SO:0000234'" in EXCLUDED
    assert "'SO:0000204'" not in EXCLUDED
    assert "'SO:0000205'" not in EXCLUDED
