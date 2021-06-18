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

import csv
import gzip
import shutil
import typing as ty
from pathlib import Path

import joblib

from rnacentral_pipeline.rnacentral.r2dt import parser
from rnacentral_pipeline.rnacentral.r2dt import should_show
from rnacentral_pipeline.rnacentral.r2dt.models import crw
from rnacentral_pipeline.rnacentral.r2dt.models import gtrnadb
from rnacentral_pipeline.rnacentral.r2dt.models import ribovision
from rnacentral_pipeline.rnacentral.r2dt.models import rnase_p


def parse(model_mapping: ty.TextIO, directory: str, allow_missing=False):
    path = Path(directory)
    return parser.parse(model_mapping, path, allow_missing=allow_missing)


def write(
    model_mapping: ty.TextIO, directory: str, output: ty.TextIO, allow_missing=False
):
    """
    Parse all the secondary structure data from the given directory and write
    it to the given file.
    """

    parsed = parse(model_mapping, directory, allow_missing=allow_missing)
    writeable = (e.writeable() for e in parsed)
    csv.writer(output).writerows(writeable)


def publish(
    model_mapping: ty.TextIO,
    directory: str,
    output: str,
    allow_missing=False,
    suffix="",
):
    out_path = Path(output)
    for result in parse(model_mapping, directory, allow_missing=allow_missing):
        publish_path = out_path / result.publish_path(suffix=suffix, compressed=True)
        try:
            publish_path.parent.mkdir(parents=True, exist_ok=True)
        except FileExistsError:
            if not publish_path.parent.exists():
                raise ValueError("Could not create publishing directory")

        with gzip.open(publish_path, "wb") as out:
            with result.info.svg.open("rb") as inp:
                shutil.copyfileobj(inp, out)


def prepare_s3(
    model_mapping: ty.TextIO,
    directory: str,
    output: Path,
    file_list: Path,
    allow_missing=False,
):

    if not output.exists():
        output.mkdir(parents=True)

    with file_list.open("w") as raw:
        seen = set()
        results = parse(model_mapping, directory, allow_missing=allow_missing)
        for result in results:
            if result.urs in seen:
                raise ValueError(f"Dupcliate URS {result.urs}")
            seen.add(result.urs)
            s3_path = output / f"{result.urs}.svg.gz"
            if s3_path.exists():
                raise ValueError(f"Will not overwrite {s3_path}")

            if not result.info.svg.exists():
                raise ValueError(f"Somehow missing unnormalized path {result.info.svg}")

            with gzip.open(s3_path, "wb") as out:
                with result.info.svg.open("rb") as inp:
                    shutil.copyfileobj(inp, out)
            raw.write(str(s3_path))
            raw.write("\n")


def write_model(generator, handle, output):
    data = generator(handle)
    data = (d.writeable() for d in data)
    csv.writer(output).writerows(data)


def write_gtrnadb(handle, output):
    return write_model(gtrnadb.parse, handle, output)


def write_ribovision(handle, output):
    return write_model(ribovision.parse, handle, output)


def write_crw(handle, output):
    return write_model(crw.parse, handle, output)


def write_rnase_p(handle, output):
    return write_model(rnase_p.parse, handle, output)


def write_should_show(model: Path, handle: ty.IO, db_url: str, output: ty.IO):
    return should_show.write(model, handle, db_url, output)


def write_training_data(handle: ty.IO, db_url: str, output: ty.IO):
    return should_show.write_training_data(handle, db_url, output)


def build_model(handle: ty.IO, db_url: str, output: Path):
    return should_show.write_model(handle, db_url, output)


def write_converted_sheet(handle: ty.IO, output: ty.IO):
    return should_show.convert_sheet(handle, output)


def write_inspect_data(handle: ty.IO, db_url: str, output: ty.IO):
    return should_show.write_inspect_data(handle, db_url, output)
