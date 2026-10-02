# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import pickle

from .. import builders as b

RRNA = "SO:0000252"
LSU_RRNA = "SO:0000651"
TRNA = "SO:0000253"


def test_so_term_for_builds_a_matching_rna_type():
    context = b.context()
    result = context.so_term_for("tRNA")
    assert result.so_term.so_id == TRNA


def test_validate_delegates_to_repeats(monkeypatch):
    context = b.context()
    called = []
    monkeypatch.setattr(context.repeats, "validate", lambda: called.append(True))
    context.validate()
    assert called == [True]


def test_dump_writes_a_loadable_pickle(tmp_path):
    context = b.context()
    out = tmp_path / "context.pickle"
    context.dump(out)

    with out.open("rb") as raw:
        data = pickle.load(raw)
    assert data == {"repeats": context.repeats}


def test_term_is_a_true_for_the_same_term():
    context = b.context()
    term = b.rna_type(RRNA)
    assert context.term_is_a("rRNA", term) is True


def test_term_is_a_true_for_an_ancestor_term():
    context = b.context()
    term = b.rna_type(LSU_RRNA)
    assert context.term_is_a("rRNA", term) is True


def test_term_is_a_false_for_an_unrelated_term():
    context = b.context()
    term = b.rna_type(TRNA)
    assert context.term_is_a("rRNA", term) is False
