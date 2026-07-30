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

"""
Deploys the in-database PL/pgSQL functions that the release pipeline relies on.

Source of truth is ``database_functions/<schema>/<func>.sql`` (one full
``CREATE OR REPLACE`` per file, as dumped by ``fetch_functions.sh``). Each file's
sha256 is recorded in a tracking table after it is applied, so re-running only
touches the functions whose source actually changed.

``CREATE OR REPLACE FUNCTION`` for PL/pgSQL is order-independent (bodies are not
resolved against other objects until run time), so files are applied in plain
(schema, name) order with no dependency graph.
"""

import dataclasses as dc
import hashlib
import logging
import typing as ty
from pathlib import Path

import psycopg2

from rnacentral_pipeline import db

LOGGER = logging.getLogger(__name__)

# database_functions/ lives at the repo root: <root>/rnacentral_pipeline/rnacentral/release/functions.py
FUNCTIONS_DIR = Path(__file__).resolve().parents[3] / "database_functions"

# Where the applied-checksum bookkeeping lives. Kept in the main rnacen schema
# alongside everything else it manages.
TRACKING_TABLE = "rnacen.applied_functions"

# Schemas that are tracked in git but must never be pushed by a release deploy --
# e.g. rnc_test holds test fixtures (setup/teardown/check_*) that have no business
# running against production.
DEPLOY_EXCLUDE_SCHEMAS = frozenset({"rnc_test"})


@dc.dataclass(frozen=True)
class Function:
    schema: str
    name: str
    path: Path
    sql: str

    @property
    def key(self) -> ty.Tuple[str, str]:
        return (self.schema, self.name)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.sql.encode("utf-8")).hexdigest()


def discover(
    base: Path = FUNCTIONS_DIR,
    exclude_schemas: ty.Collection[str] = (),
) -> ty.List[Function]:
    """
    Load every <schema>/<func>.sql file under base, sorted by (schema, name).

    Schemas in exclude_schemas are skipped -- used to keep deploy-excluded schemas
    (e.g. rnc_test) tracked in git without applying them to the database.
    """
    exclude = set(exclude_schemas)
    functions = []
    for path in sorted(base.glob("*/*.sql")):
        if path.parent.name in exclude:
            continue
        functions.append(
            Function(
                schema=path.parent.name,
                name=path.stem,
                path=path,
                sql=path.read_text(),
            )
        )
    if not functions:
        raise ValueError(f"No function files found under {base}")
    return functions


def _applied_shas(cur) -> ty.Dict[ty.Tuple[str, str], str]:
    cur.execute(f"SELECT schema_name, function_name, sha256 FROM {TRACKING_TABLE}")
    return {(s, n): sha for (s, n, sha) in cur.fetchall()}


def pending(cur, functions: ty.List[Function]) -> ty.List[Function]:
    """Functions whose source differs from what was last applied (new or changed)."""
    applied = _applied_shas(cur)
    return [f for f in functions if applied.get(f.key) != f.sha256]


def apply(db_url: str, dry_run: bool = False) -> ty.List[Function]:
    """
    Apply every function whose source changed since it was last applied.

    All changes run in a single transaction: either every pending function is
    replaced and recorded, or nothing is. Returns the list that was (or would be)
    applied.
    """
    functions = discover(exclude_schemas=DEPLOY_EXCLUDE_SCHEMAS)
    with db.connect(db_url) as conn:
        with conn.cursor() as cur:
            todo = pending(cur, functions)
            if not todo:
                LOGGER.info("All %d functions up to date", len(functions))
                return []

            for fn in todo:
                LOGGER.info(
                    "%s %s.%s",
                    "Would apply" if dry_run else "Applying",
                    fn.schema,
                    fn.name,
                )
                if dry_run:
                    continue
                cur.execute(fn.sql)
                cur.execute(
                    f"""
                    INSERT INTO {TRACKING_TABLE}
                        (schema_name, function_name, sha256, applied_at)
                    VALUES (%s, %s, %s, clock_timestamp())
                    ON CONFLICT (schema_name, function_name)
                    DO UPDATE SET sha256 = excluded.sha256,
                                  applied_at = excluded.applied_at
                    """,
                    (fn.schema, fn.name, fn.sha256),
                )

        if dry_run:
            conn.rollback()
        else:
            conn.commit()

    LOGGER.info("%s %d function(s)", "Would apply" if dry_run else "Applied", len(todo))
    return todo
