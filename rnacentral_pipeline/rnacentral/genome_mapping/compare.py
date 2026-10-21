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

import gzip
import logging
import re
import typing as ty
from pathlib import Path

import polars as pl
import psycopg2 as pg

LOGGER = logging.getLogger(__name__)

URS_REGEX = re.compile(r"(URS[0-9A-F]{10}_\d+)")

GLOBS = {"gff": "*.gff3.gz", "bed": "*.bed.gz"}


def written_ids(ftp_path: ty.Union[str, Path], kind: str) -> pl.DataFrame:
    """
    Collect every URS_taxid written to the exported files of the given kind
    ('gff' or 'bed') in ftp_path.
    """
    all_ids: ty.Set[str] = set()
    for path in sorted(Path(ftp_path).glob(GLOBS[kind])):
        file_ids: ty.Set[str] = set()
        coordinates = 0
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                match = URS_REGEX.search(line)
                if match:
                    coordinates += 1
                    file_ids.add(match.group(1))
        LOGGER.info(
            "%s: %i coordinates, %i unique ids", path.name, coordinates, len(file_ids)
        )
        all_ids.update(file_ids)
    LOGGER.info("Total written ids: %i", len(all_ids))
    return pl.DataFrame({"urs_taxid": list(all_ids)}, schema={"urs_taxid": pl.String})


def database_ids(db_url: str) -> pl.DataFrame:
    """
    Get every active URS_taxid that the database says should have coordinates.
    """
    query = "SELECT id FROM rnc_rna_precomputed WHERE is_active AND has_coordinates"
    conn = pg.connect(db_url)
    try:
        with conn.cursor() as cur:
            cur.execute(query)
            rows = cur.fetchall()
    finally:
        conn.close()
    return pl.DataFrame(
        {"urs_taxid": [r[0] for r in rows]}, schema={"urs_taxid": pl.String}
    ).unique()


def missing_from_files(
    ftp_path: ty.Union[str, Path], db_url: str, kind: str
) -> pl.DataFrame:
    """
    Find the URS_taxids the database expects to have coordinates that are
    absent from the exported files.
    """
    in_files = written_ids(ftp_path, kind)
    in_db = database_ids(db_url)
    missing = in_db.join(in_files, on="urs_taxid", how="anti")
    LOGGER.info(
        "%i ids in files, %i in database, %i in database but not in files",
        len(in_files),
        len(in_db),
        len(missing),
    )
    return missing
