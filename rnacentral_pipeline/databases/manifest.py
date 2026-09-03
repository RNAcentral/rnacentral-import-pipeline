# -*- coding: utf-8 -*-

"""
Shared machinery for incremental ("delta") parsing.

A database's import can skip the expensive per-record parsing for records that
have not changed since the last import. To decide what changed we keep a manifest
of ``accession -> signature`` (a hash of the raw record) from the last successful
import, diff the new input against it, and only fully parse the new/changed
records. Records that vanished from the input become an explicit deletion list.

This module is database-agnostic: signature/diff are pure, and the manifest table is
a generic ``pipeline_tracking_import`` keyed by (database, accession), LIST-partitioned
on ``database`` so ENA's bulk never drives another database's vacuum or bloat.
See docs/incremental-parsing.md.
"""

import csv
import hashlib
import json
import typing as ty
from pathlib import Path

import attr
import psycopg2
from psycopg2 import sql
from psycopg2.extras import execute_values

MANIFEST_TABLE = "rnacen.pipeline_tracking_import"

MANIFEST_CSV = "manifest.csv"
DELETIONS_CSV = "deletions.csv"

# database leads the PK because every query scopes to one, and a partition key must
# be part of it.
CREATE_MANIFEST_SQL = f"""
CREATE TABLE IF NOT EXISTS {MANIFEST_TABLE} (
    database   text        NOT NULL,
    accession  text        NOT NULL,
    signature  text        NOT NULL,
    updated_at timestamptz NOT NULL DEFAULT clock_timestamp(),
    PRIMARY KEY (database, accession)
) PARTITION BY LIST (database)
"""


def partition_name(database: str) -> str:
    """
    Partition identifier for a database, e.g. ENA -> pipeline_tracking_import_ena.
    Normalised here because Postgres lower-cases identifiers and caps them at 63 bytes.
    """
    safe = "".join(c if c.isalnum() else "_" for c in database.lower())
    return f"pipeline_tracking_import_{safe}"[:63]


def record_signature(raw: ty.Any) -> str:
    """
    A stable hash of a raw record. Canonical JSON (sorted keys) so the signature
    depends only on content, not key order or whitespace. The whole record is
    hashed, so any change to any field is caught.
    """
    canonical = json.dumps(raw, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def text_signature(raw: str) -> str:
    """
    A stable hash of an already-serialised raw record (e.g. an EMBL record block).
    Used where the raw form is text rather than a structured object; hashing the
    whole text can never produce a false "unchanged".
    """
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


@attr.s(frozen=True)
class Diff:
    """The four disjoint outcomes of comparing a new manifest against the old one."""

    new = attr.ib(type=ty.FrozenSet[str])
    changed = attr.ib(type=ty.FrozenSet[str])
    dropped = attr.ib(type=ty.FrozenSet[str])
    unchanged = attr.ib(type=ty.FrozenSet[str])

    @property
    def to_parse(self) -> ty.FrozenSet[str]:
        """Accessions that must be fully parsed this import (new or changed)."""
        return self.new | self.changed

    @property
    def is_bootstrap(self) -> bool:
        """True when there was no prior manifest (everything is new)."""
        return not (self.changed or self.dropped or self.unchanged)


def compute_diff(
    new_signatures: ty.Mapping[str, str],
    old_signatures: ty.Mapping[str, str],
) -> Diff:
    """Partition accessions into new / changed / dropped / unchanged."""
    new_keys = set(new_signatures)
    old_keys = set(old_signatures)

    new = new_keys - old_keys
    dropped = old_keys - new_keys
    common = new_keys & old_keys
    changed = {k for k in common if new_signatures[k] != old_signatures[k]}
    unchanged = common - changed

    return Diff(
        new=frozenset(new),
        changed=frozenset(changed),
        dropped=frozenset(dropped),
        unchanged=frozenset(unchanged),
    )


def ensure_table(conn) -> None:
    """Create the partitioned parent table if it does not exist yet."""
    with conn.cursor() as cur:
        cur.execute(CREATE_MANIFEST_SQL)
    conn.commit()


def ensure_partition(conn, database: str) -> None:
    """
    Create this database's partition if missing. Required before any write; reads are
    fine without one, matching nothing, which is the right answer for a new database.
    """
    ensure_table(conn)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                "CREATE TABLE IF NOT EXISTS rnacen.{} PARTITION OF {} FOR VALUES IN (%s)"
            ).format(
                sql.Identifier(partition_name(database)),
                sql.SQL(MANIFEST_TABLE),
            ),
            (database,),
        )
    conn.commit()


def load_signatures(conn, database: str) -> ty.Dict[str, str]:
    """The stored ``accession -> signature`` map for a database (empty if none)."""
    ensure_table(conn)
    with conn.cursor() as cur:
        cur.execute(
            f"SELECT accession, signature FROM {MANIFEST_TABLE} WHERE database = %s",
            (database,),
        )
        return dict(cur.fetchall())


def load_signatures_for(db_url: str, database: str) -> ty.Dict[str, str]:
    """Convenience wrapper that opens its own connection."""
    conn = psycopg2.connect(db_url)
    try:
        return load_signatures(conn, database)
    finally:
        conn.close()


def store_signatures(
    conn,
    database: str,
    signatures: ty.Mapping[str, str],
    dropped: ty.Iterable[str] = (),
) -> None:
    """
    Replace the stored manifest for a database: upsert every current signature and
    remove dropped accessions. Call this only after the database's load succeeds,
    so a failed load leaves the previous manifest intact.
    """
    ensure_partition(conn, database)
    dropped = list(dropped)
    with conn.cursor() as cur:
        if dropped:
            cur.execute(
                f"DELETE FROM {MANIFEST_TABLE} WHERE database = %s AND accession = ANY(%s)",
                (database, dropped),
            )
        # Batched upsert: a per-row INSERT loop is fine for a few thousand HGNC rows
        # but hopeless for ENA's millions, so send them in pages.
        rows = (
            (database, accession, signature)
            for accession, signature in signatures.items()
        )
        # Without the WHERE, every run rewrites every row -- at ENA scale hundreds of
        # millions of dead tuples for a handful of changed signatures.
        execute_values(
            cur,
            f"""
            INSERT INTO {MANIFEST_TABLE} AS m (database, accession, signature)
            VALUES %s
            ON CONFLICT (database, accession)
            DO UPDATE SET signature = excluded.signature,
                          updated_at = clock_timestamp()
            WHERE m.signature IS DISTINCT FROM excluded.signature
            """,
            rows,
            page_size=5000,
        )
    conn.commit()


def write_artifacts(
    output_dir: ty.Union[str, Path],
    database: str,
    signatures: ty.Mapping[str, str],
    deletions: ty.Iterable[str],
) -> None:
    """
    Write the two side-channel files a delta parse produces, for the load step:
      * manifest.csv  -- database,accession,signature for every current record;
      * deletions.csv -- database,accession for every dropped record.
    """
    out = Path(output_dir)
    with (out / MANIFEST_CSV).open("w", newline="") as handle:
        writer = csv.writer(handle)
        for accession, signature in sorted(signatures.items()):
            writer.writerow([database, accession, signature])
    with (out / DELETIONS_CSV).open("w", newline="") as handle:
        writer = csv.writer(handle)
        for accession in deletions:
            writer.writerow([database, accession])


def apply_artifacts(
    conn,
    manifest_csv: ty.Union[str, Path],
    deletions_csv: ty.Optional[ty.Union[str, Path]] = None,
) -> None:
    """
    Promote a delta parse's manifest into pipeline_tracking_import. Call only after the
    database's load/release has committed. Reads the (database, accession, signature)
    rows from manifest_csv and the dropped (database, accession) rows from
    deletions_csv, then replaces each database's stored signatures.
    """
    signatures: ty.Dict[str, ty.Dict[str, str]] = {}
    with Path(manifest_csv).open("r", newline="") as handle:
        for database, accession, signature in csv.reader(handle):
            signatures.setdefault(database, {})[accession] = signature

    dropped: ty.Dict[str, ty.List[str]] = {}
    if deletions_csv is not None and Path(deletions_csv).exists():
        with Path(deletions_csv).open("r", newline="") as handle:
            for database, accession in csv.reader(handle):
                dropped.setdefault(database, []).append(accession)

    for database in set(signatures) | set(dropped):
        store_signatures(
            conn,
            database,
            signatures.get(database, {}),
            dropped.get(database, []),
        )


class _CopySource:
    """
    A read-only file object that serialises ``(accession, signature)`` tuples to CSV
    lazily, one row at a time, for ``copy_expert``. At ENA scale the full signature
    set is multiple GB, so buffering it into a single StringIO before COPY is what
    OOM-killed the diff; this pulls rows from the iterator only as psycopg2 reads.
    """

    def __init__(self, rows: ty.Iterable[ty.Tuple[str, str]]):
        self._writer = csv.writer(self)
        self._rows = iter(rows)
        self._buf = ""
        self._out = ""

    def write(self, text: str) -> int:
        # csv.writer writes here; stash the rendered row for read() to hand out.
        self._out += text
        return len(text)

    def read(self, size: int = -1) -> str:
        while size < 0 or len(self._buf) < size:
            try:
                self._writer.writerow(next(self._rows))
            except StopIteration:
                break
            self._buf += self._out
            self._out = ""
        if size < 0:
            chunk, self._buf = self._buf, ""
            return chunk
        chunk, self._buf = self._buf[:size], self._buf[size:]
        return chunk


@attr.s(frozen=True)
class DbDiff:
    """Result of a database-side diff (see :func:`diff_via_db`)."""

    to_parse = attr.ib(type=ty.List[str])  # accessions new or changed since last import
    deletions = attr.ib(type=ty.List[str])  # accessions present last import, now absent
    is_bootstrap = attr.ib(type=bool)  # True when there was no prior manifest


def diff_via_db(
    conn,
    database: str,
    signatures: ty.Iterable[ty.Tuple[str, str]],
) -> DbDiff:
    """
    Diff a large new signature set against the stored manifest inside Postgres.

    ``signatures`` is an iterable of ``(accession, signature)`` for every record in
    the new import. For databases the size of ENA, loading both the new and the old
    manifest into Python to diff them is memory-hungry; instead the new set is
    COPYed into a temp table and the set differences are computed in SQL against
    :data:`MANIFEST_TABLE` scoped to ``database``.

    Returns the accessions to parse (new or changed), the accessions to delete
    (stored but now absent), and whether this is a bootstrap (no prior manifest for
    the database). On bootstrap ``to_parse`` is empty and callers should parse
    everything -- the caller, not this function, decides how to represent "all".
    """
    ensure_table(conn)
    with conn.cursor() as cur:
        # Stream into an unconstrained staging table first: ENA can emit the same
        # location-based accession twice within one snapshot, which would abort a
        # COPY straight into a PRIMARY KEY table. Collapse the duplicates on the way
        # in, the same last-wins way store_signatures/apply_artifacts already do.
        cur.execute(
            "CREATE TEMP TABLE _new_manifest_raw (accession text NOT NULL, signature text NOT NULL) ON COMMIT DROP"
        )

        cur.copy_expert(
            "COPY _new_manifest_raw (accession, signature) FROM STDIN WITH CSV",
            _CopySource(signatures),
        )

        cur.execute(
            "CREATE TEMP TABLE _new_manifest (accession text PRIMARY KEY, signature text NOT NULL) ON COMMIT DROP"
        )
        cur.execute(
            "INSERT INTO _new_manifest (accession, signature) "
            "SELECT DISTINCT ON (accession) accession, signature "
            "FROM _new_manifest_raw ORDER BY accession"
        )

        cur.execute(
            f"SELECT count(*) FROM {MANIFEST_TABLE} WHERE database = %s", (database,)
        )
        is_bootstrap = cur.fetchone()[0] == 0

        to_parse: ty.List[str] = []
        if not is_bootstrap:
            cur.execute(
                f"""
                SELECT n.accession
                FROM _new_manifest n
                LEFT JOIN {MANIFEST_TABLE} m
                  ON m.database = %s AND m.accession = n.accession
                WHERE m.accession IS NULL OR m.signature <> n.signature
                """,
                (database,),
            )
            to_parse = [row[0] for row in cur]

        cur.execute(
            f"""
            SELECT m.accession
            FROM {MANIFEST_TABLE} m
            LEFT JOIN _new_manifest n ON n.accession = m.accession
            WHERE m.database = %s AND n.accession IS NULL
            """,
            (database,),
        )
        deletions = [row[0] for row in cur]

    return DbDiff(to_parse=to_parse, deletions=deletions, is_bootstrap=is_bootstrap)
