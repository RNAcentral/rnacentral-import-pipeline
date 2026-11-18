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
"""

import csv

import pyarrow.parquet as pq
import pytest

from rnacentral_pipeline import schemas, writers
from rnacentral_pipeline.databases import data

COLUMNS = [
    "accession",
    "assigned_by",
    "disease_text",
    "ontology_term_id",
    "match_type",
    "pubmed_id",
]

MAPPED = ["ACC.1", "EVLNCRNAS", "gallbladder cancer", "MONDO:0005411", "label_exact"]
UNMAPPED = ["ACC.1", "EVLNCRNAS", "GBC", None, "unmapped"]


def entry():
    return data.Entry(
        primary_id="primary",
        accession="ACC.1",
        ncbi_tax_id=9606,
        database="EVLNCRNAS",
        sequence="ACGTACGTACGTACGT",
        regions=[],
        rna_type="SO:0001877",
        url="https://example.org/ACC.1",
        seq_version="1",
        description="A lncRNA",
        disease_annotations=[
            data.DiseaseAnnotation(
                disease_text="gallbladder cancer",
                ontology_term_id="MONDO:0005411",
                match_type="label_exact",
                pubmed_ids=[31174563, 123],
                assigned_by="EVLNCRNAS",
            ),
            data.DiseaseAnnotation(
                disease_text="GBC",
                ontology_term_id=None,
                match_type="unmapped",
                pubmed_ids=[],
                assigned_by="EVLNCRNAS",
            ),
        ],
    )


def test_schema_has_the_contract_column_order():
    assert schemas.DISEASE_ANNOTATIONS.names == COLUMNS


def test_csv_writer_writes_annotations_and_terms(tmp_path):
    with writers.build(writers.EntryWriter, tmp_path) as writer:
        writer.write([entry()])

    with (tmp_path / "disease_annotations.csv").open() as handle:
        rows = list(csv.reader(handle))
    # csv has no null; pgloader reads the empty field as NULL.
    assert rows == [
        MAPPED + ["31174563"],
        MAPPED + ["123"],
        ["" if v is None else v for v in UNMAPPED] + [""],
    ]

    with (tmp_path / "terms.csv").open() as handle:
        terms = [row[0] for row in csv.reader(handle)]
    assert terms == ["SO:0001877", "MONDO:0005411"]


def test_parquet_writer_writes_annotations_and_terms(tmp_path):
    with writers.parquet_entry_writer(tmp_path) as writer:
        writer.write([entry()])

    table = pq.read_table(tmp_path / "disease_annotations.parquet")
    assert table.schema.equals(schemas.DISEASE_ANNOTATIONS)
    assert [list(row.values()) for row in table.to_pylist()] == [
        MAPPED + ["31174563"],
        MAPPED + ["123"],
        UNMAPPED + [None],
    ]

    terms = pq.read_table(tmp_path / "terms.parquet").column(0).to_pylist()
    assert terms == ["SO:0001877", "MONDO:0005411"]


@pytest.mark.parametrize(
    "term,match_type",
    [
        (None, "label_exact"),
        ("MONDO:0005411", "unmapped"),
        ("MONDO:0005411", "fuzzy"),
    ],
)
def test_rejects_inconsistent_term_and_match_type(term, match_type):
    with pytest.raises(ValueError):
        data.DiseaseAnnotation(
            disease_text="gallbladder cancer",
            ontology_term_id=term,
            match_type=match_type,
            pubmed_ids=[],
            assigned_by="EVLNCRNAS",
        )
