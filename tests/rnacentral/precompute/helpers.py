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

import os
import typing as ty

from rnacentral_pipeline.databases.sequence_ontology import tree
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from tests.helpers import run_range_as_single, run_with_replacements

SO_TREE = tree.load_ontology(tree.REMOTE_ONTOLOGY)


def load_data(rna_id: str) -> ty.Tuple[Context, Sequence]:
    path = os.path.join("files", "precompute", "query.sql")
    upi, taxid = rna_id.split("_")
    data = run_with_replacements(
        path,
        (":tablename", "rna"),
        (
            "todo.id BETWEEN :min AND :max",
            "xref.upi ='%s' AND xref.taxid = %i" % (upi, int(taxid)),
        ),
    )
    return (Context(so_tree=SO_TREE), Sequence.build(SO_TREE, data))


def load_for_upi(upi: str) -> ty.List[Sequence]:
    path = os.path.join("files", "precompute", "query.sql")
    loaded = list(
        run_with_replacements(
            path,
            (":tablename", "rna"),
            ("todo.id BETWEEN :min AND :max", "xref.upi ='%s'" % upi),
            take_all=True,
        )
    )
    return [Sequence.build(SO_TREE, l) for l in loaded]


def load_for_range(start: int, stop: int) -> ty.List[Sequence]:
    path = os.path.join("files", "precompute", "query.sql")
    loaded = run_with_replacements(
        path,
        (":min", str(start)),
        (":max", str(stop)),
        (":tablename", "rna"),
        take_all=True,
    )
    return [Sequence.build(SO_TREE, l) for l in loaded]
