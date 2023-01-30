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
"""

import csv
import enum
import logging
import typing as ty
from pathlib import Path

import joblib
import pandas as pd
import psycopg2
import psycopg2.extras
from more_itertools import chunked
from pypika import Query, Table
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, train_test_split

LOGGER = logging.getLogger(__name__)


SOURCE_MAP = {
    "crw": 0,
    "ribovision": 1,
    "gtrnadb": 2,
    "rnase_p": 3,
    "rfam": 4,
}


@enum.unique
class Attributes(enum.Enum):
    SourceIndex = "source_index"
    SequenceLength = "sequence_length"
    DiagramSequenceLength = "diagram_sequence_length"
    ModelLength = "model_length"
    ModelBasepairCount = "model_basepair_count"
    DiagramBps = "diagram_bps"
    DiagramModelLength = "diagram_model_length"
    DiagramOverlapCount = "diagram_overlap_count"

    @classmethod
    def model_columns(cls) -> ty.List[str]:
        return [attr.column_name() for attr in cls]

    def column_name(self) -> str:
        return self.value


MODEL_COLUMNS: ty.List[str] = Attributes.model_columns()


def chunked_query(
    ids: ty.Iterable[str], query_builder, db_url: str, chunk_size=100
) -> ty.Iterable[ty.Dict[str, ty.Any]]:
    conn = psycopg2.connect(db_url)
    for chunk in chunked(ids, chunk_size):
        sql = str(query_builder(chunk))
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
            cur.execute(sql)
            for result in cur:
                yield dict(result)


def fetch_modeled_data(
    all_ids: ty.Iterable[str], db_url: str, chunk_size=100
) -> ty.Iterable[ty.Dict[str, ty.Any]]:
    rna = Table("rna")
    ss = Table("rnc_secondary_structure_layout")
    sm = Table("rnc_secondary_structure_layout_models")

    def build_query(ids):
        return (
            Query.from_(rna)
            .select(
                rna.upi.as_("urs"),
                rna.len.as_("sequence_length"),
                sm.model_source,
                ss.sequence_start.as_("diagram_sequence_start"),
                ss.sequence_stop.as_("diagram_sequence_stop"),
                ss.basepair_count.as_("diagram_bps"),
                ss.model_start.as_("diagram_model_start"),
                ss.model_stop.as_("diagram_model_stop"),
                sm.model_length,
                sm.model_basepair_count,
                ss.overlap_count.as_("diagram_overlap_count"),
            )
            .join(ss)
            .on(ss.urs == rna.upi)
            .join(sm)
            .on(sm.id == ss.model_id)
            .where(ss.urs.isin(ids))
        )

    seen: ty.Set[str] = set()
    results = chunked_query(all_ids, build_query, db_url, chunk_size=chunk_size)
    for result in results:
        if any(v is None for v in result.values()):
            continue
        yield result
        seen.add(result["urs"])

    for urs in all_ids:
        if urs not in seen:
            LOGGER.warn("Missed loading %s", urs)


def infer_columns(frame: pd.DataFrame):
    frame["diagram_sequence_length"] = (
        frame["diagram_sequence_stop"] - frame["diagram_sequence_start"]
    )
    frame["diagram_model_length"] = (
        frame["diagram_model_stop"] - frame["diagram_model_start"]
    )
    frame["source_index"] = frame.model_source.map(SOURCE_MAP)
    if frame["source_index"].isnull().any():
        raise ValueError("Could not build source_index for all training data")


def fetch_training_data(handle: ty.IO, db_url: str) -> pd.DataFrame:
    ids = []
    training = {}
    for (urs, flag) in csv.reader(handle):
        ids.append(urs)
        if flag == "1":
            training[urs] = True
        elif flag == "0":
            training[urs] = False
        else:
            raise ValueError(f"Unknown flag {flag}")

    filled = []
    for metadata in fetch_modeled_data(ids, db_url):
        urs = metadata["urs"]
        if urs not in training:
            raise ValueError(f"Got an extra entry, somehow {metadata}")
        metadata["valid"] = training[urs]
        filled.append(metadata)

    training = pd.DataFrame.from_records(filled)
    infer_columns(training)
    return training


def train(handle, db_url, cross_validation=5, test_size=0.4) -> RandomForestClassifier:
    data = fetch_training_data(handle, db_url)
    X_train, X_test, y_train, y_test = train_test_split(
        data[MODEL_COLUMNS].to_numpy(), data["valid"].to_numpy(), test_size=test_size
    )

    clf = RandomForestClassifier(min_samples_split=5)
    scores = cross_val_score(clf, X_train, y_train, cv=cross_validation)
    LOGGER.info("%s fold cross validation scores: %s", cross_validation, scores)
    clf.fit(X_train, y_train)
    LOGGER.info("Test data (%f) scoring %s", test_size, clf.score(X_test, y_test))
    return clf


def from_result(clf, result) -> bool:
    predictable = {}
    for attribute in Attributes:
        value = attribute.r2dt_result_value(result)
        predictable[attribute.column_name()] = [value]
    predictable = pd.DataFrame.from_records(predictable)
    return clf.predict(predictable)[0]


def write(model_path: Path, handle: ty.IO, db_url: str, output: ty.IO):
    model = joblib.load(model_path)
    ids = [r[0] for r in csv.reader(handle)]
    modeled = fetch_modeled_data(ids, db_url)
    frame = pd.DataFrame.from_records(modeled)
    if len(frame) > 0:
        infer_columns(frame)
        predicted = model.predict(frame[MODEL_COLUMNS].to_numpy())
        to_write = pd.DataFrame()
        to_write["urs"] = frame["urs"]
        to_write["should_show"] = predicted.astype(int)
        to_write.to_csv(output, index=False)


def write_model(handle: ty.IO, db_url: str, output: Path):
    joblib.dump(train(handle, db_url), output)


def write_training_data(handle: ty.IO, db_url: str, output: ty.IO):
    ids = []
    for row in csv.reader(handle):
        ids.append(row[0])
    modeled = list(fetch_modeled_data(ids, db_url))
    writer = csv.DictWriter(output, fieldnames=modeled[0].keys())
    writer.writeheader()
    writer.writerows(modeled)


def convert_sheet(handle: ty.IO, output: ty.IO):
    converted = []
    for row in csv.DictReader(handle):
        urs = row["urs"]
        raw_should_show = row["Labeled Should show"]
        if not raw_should_show:
            LOGGER.info("No value for %s", urs)

        should_show = None
        raw_should_show = raw_should_show.lower()
        if raw_should_show == "true":
            should_show = "1"
        elif raw_should_show == "false":
            should_show = "0"
        else:
            LOGGER.warn("Unknown should show in %s", row)
            continue
        converted.append((urs, should_show))
    converted.sort(key=lambda r: r[0])
    writer = csv.writer(output)
    writer.writerows(converted)


def inspect_data(data, db_url: str) -> ty.Iterable[ty.Dict[str, ty.Any]]:
    def build_query(ids):
        ss = Table("rnc_secondary_structure_layout")
        sm = Table("rnc_secondary_structure_layout_models")
        pre = Table("rnc_rna_precomputed")
        return (
            Query.from_(ss)
            .join(sm)
            .on(sm.id == ss.model_id)
            .join(pre)
            .on(pre.urs == sm.urs)
            .select(
                sm.model_source,
                sm.model_name,
                sm.model_so_term,
            )
            .where(ss.urs.isin(ids))
            .where(pre.taxid.isnotnull)
        )

    mapping = {d[0]: d for d in data}
    seen: ty.Set[str] = set()
    results = chunked_query(data, build_query, db_url)
    for result in results:
        if any(v is None for v in result.values()):
            continue
        yield {
            "urs": result["urs"],
            "link": f"https://rnacentral.org/rna/{result['urs']}",
            "model_source": result["model_source"],
            "model_name": result["model_name"],
            "model_so_term": result["model_so_term"],
            "Labeled Should show": result["urs"],
        }
        seen.add(result["urs"])

    for urs in mapping.keys():
        if urs not in seen:
            LOGGER.warn("Missed loading %s", urs)


def write_inspect_data(handle: ty.IO, db_url: str, output: ty.IO):
    data = list(csv.reader(handle))
    inspect = list(inspect_data(data, db_url))
    writer = csv.DictWriter(output, fieldnames=inspect[0].keys())
    writer.writeheader()
    writer.writerows(inspect)
