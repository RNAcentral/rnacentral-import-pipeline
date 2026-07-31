# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

The classifier tests run against data/r2dt/should-show/labelled-corpus.csv, a
frozen snapshot of should_show.fetch_modeled_data output for the labelled URS in
training-data.csv, plus the label and the shipped model's prediction. It is
snapshotted rather than queried so the tests are deterministic and need no
database: the pipeline rewrites r2dt_results on every re-scan, so a live query
would make these results depend on when they were run.

Two URS in training-data.csv carry contradictory labels (URS000000BB03 and
URS000028467F appear as both 0 and 1) and are excluded from the corpus.

The committed corpus was generated from pfmrnac16pro. To regenerate after
retraining, run fetch_modeled_data over the training URS against a populated
database, add the label and prediction columns, and update the constants below
with the new measurements -- and note which database it came from, since the
feature values differ between instances.
"""

import csv
import io
from pathlib import Path

import joblib
import pandas as pd
import pytest

from rnacentral_pipeline.rnacentral.r2dt import should_show

MODEL = Path("files/r2dt/should-show/model.joblib")
CORPUS = Path("data/r2dt/should-show/labelled-corpus.csv")

# Measured against the shipped model; a change here is a change in behaviour,
# not a flaky threshold. See the module docstring for how to regenerate.
EXPECTED_ACCURACY = 0.8965
WRONGLY_SHOWN = 43
WRONGLY_HIDDEN = 148


@pytest.mark.r2dt
def test_model_columns_are_the_attribute_values():
    assert should_show.Attributes.model_columns() == [
        "source_index",
        "sequence_length",
        "diagram_sequence_length",
        "model_length",
        "model_basepair_count",
        "diagram_bps",
        "diagram_model_length",
        "diagram_overlap_count",
    ]


@pytest.mark.r2dt
def test_infer_columns_derives_lengths_and_source_index():
    frame = pd.DataFrame.from_records(
        [
            {
                "model_source": "crw",
                "diagram_sequence_start": 3,
                "diagram_sequence_stop": 1588,
                "diagram_model_start": 3,
                "diagram_model_stop": 1512,
            },
            {
                "model_source": "rfam",
                "diagram_sequence_start": 1,
                "diagram_sequence_stop": 95,
                "diagram_model_start": 1,
                "diagram_model_stop": 108,
            },
        ]
    )
    should_show.infer_columns(frame)
    assert list(frame["diagram_sequence_length"]) == [1585, 94]
    assert list(frame["diagram_model_length"]) == [1509, 107]
    assert list(frame["source_index"]) == [0, 4]


@pytest.mark.r2dt
def test_infer_columns_rejects_an_unknown_source():
    frame = pd.DataFrame.from_records(
        [
            {
                "model_source": "not-a-source",
                "diagram_sequence_start": 1,
                "diagram_sequence_stop": 2,
                "diagram_model_start": 1,
                "diagram_model_stop": 2,
            }
        ]
    )
    with pytest.raises(ValueError):
        should_show.infer_columns(frame)


@pytest.mark.r2dt
def test_convert_sheet_maps_labels_and_sorts():
    raw = io.StringIO(
        "urs,Labeled Should show\n"
        "URS0000000002,TRUE\n"
        "URS0000000001,false\n"
        "URS0000000003,maybe\n"
        "URS0000000004,\n"
    )
    out = io.StringIO()
    should_show.convert_sheet(raw, out)
    assert out.getvalue().splitlines() == [
        "URS0000000001,0",
        "URS0000000002,1",
    ]


def corpus():
    with CORPUS.open("r") as raw:
        rows = list(csv.DictReader(raw))
    ints = {
        "sequence_length",
        "diagram_sequence_start",
        "diagram_sequence_stop",
        "diagram_bps",
        "diagram_model_start",
        "diagram_model_stop",
        "model_length",
        "model_basepair_count",
        "diagram_overlap_count",
    }
    records = [
        {
            k: (int(v) if k in ints else v)
            for k, v in r.items()
            if k in ints or k == "urs" or k == "model_source"
        }
        for r in rows
    ]
    labels = [int(r["label"]) for r in rows]
    predictions = [int(r["prediction"]) for r in rows]
    return records, labels, predictions


@pytest.mark.r2dt
def test_the_pickled_model_still_loads():
    # The pickle was written under sklearn 0.24.1 and does not load on >=1.3, so
    # a dependency bump silently breaks `should-show compute` at joblib.load.
    model = joblib.load(MODEL)
    assert model.n_features_in_ == len(should_show.MODEL_COLUMNS)


@pytest.mark.r2dt
def test_model_predictions_are_unchanged(monkeypatch):
    """Every prediction over the labelled corpus, through the production path.

    A change detector, not a correctness claim: if a sklearn upgrade or a
    retrained model moves any prediction, the diff is what needs reviewing.
    """
    records, _, predictions = corpus()
    monkeypatch.setattr(should_show, "fetch_modeled_data", lambda *a, **k: records)

    urs_list = io.StringIO("\n".join(r["urs"] for r in records) + "\n")
    out = io.StringIO()
    should_show.write(MODEL, urs_list, "postgres://unused", out)

    written = [(row[0], int(row[1])) for row in csv.reader(io.StringIO(out.getvalue()))]
    assert written == [(r["urs"], p) for r, p in zip(records, predictions)]


@pytest.mark.r2dt
def test_model_accuracy_against_the_labels():
    """Guards the quality of the shipped model, not just its stability.

    The heuristic in data.ShowInfo.showable() scores 0.867 on the same corpus,
    with 140 wrongly shown; the model is both more accurate and biased towards
    hiding, which is the safer error for a public site.
    """
    records, labels, _ = corpus()
    frame = pd.DataFrame.from_records(records)
    should_show.infer_columns(frame)
    model = joblib.load(MODEL)
    predicted = model.predict(frame[should_show.MODEL_COLUMNS].to_numpy()).astype(int)

    correct = sum(p == l for p, l in zip(predicted, labels))
    assert correct / len(labels) == pytest.approx(EXPECTED_ACCURACY, abs=0.0001)
    assert sum(1 for p, l in zip(predicted, labels) if p and not l) == WRONGLY_SHOWN
    assert sum(1 for p, l in zip(predicted, labels) if not p and l) == WRONGLY_HIDDEN


@pytest.mark.r2dt
def test_write_produces_nothing_when_there_are_no_rows(monkeypatch):
    monkeypatch.setattr(should_show, "fetch_modeled_data", lambda *a, **k: [])
    out = io.StringIO()
    should_show.write(MODEL, io.StringIO(""), "postgres://unused", out)
    assert out.getvalue() == ""


@pytest.mark.r2dt
def test_corpus_labels_are_consistent_with_training_data():
    """The corpus must stay a faithful subset of the committed training data."""
    seen = {}
    with open("data/r2dt/should-show/training-data.csv", "r") as raw:
        for urs, flag in csv.reader(raw):
            seen.setdefault(urs, set()).add(flag)

    rows = list(csv.DictReader(CORPUS.open("r")))
    assert {r["urs"] for r in rows} <= set(seen)
    for row in rows:
        assert seen[row["urs"]] == {row["label"]}, row["urs"]
