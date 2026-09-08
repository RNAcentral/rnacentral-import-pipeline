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
import polars as pl
import psycopg2
from psycopg2 import sql
from psycopg2.extras import execute_values

MANIFEST_TABLE = "rnacen.pipeline_tracking_import"
FILES_TABLE = "rnacen.pipeline_tracking_import_files"

MANIFEST_CSV = "manifest.csv"
DELETIONS_CSV = "deletions.csv"
SOURCES_CSV = "sources.csv"

# database leads the PK because every query scopes to one, and a partition key must
# be part of it.
CREATE_MANIFEST_SQL = f"""
CREATE TABLE IF NOT EXISTS {MANIFEST_TABLE} (
    database   text        NOT NULL,
    accession  text        NOT NULL,
    signature  text        NOT NULL,
    updated_at timestamptz NOT NULL DEFAULT clock_timestamp(),
    source_id  bigint,
    PRIMARY KEY (database, accession)
) PARTITION BY LIST (database)
"""

# A surrogate id rather than this convention's usual natural key: the manifest holds
# hundreds of millions of ENA rows, and 8 bytes per row beats repeating a ~60 byte
# path in every one of them.
CREATE_FILES_SQL = f"""
CREATE TABLE IF NOT EXISTS {FILES_TABLE} (
    id         bigserial   PRIMARY KEY,
    database   text        NOT NULL,
    path       text        NOT NULL,
    signature  text        NOT NULL,
    updated_at timestamptz NOT NULL DEFAULT clock_timestamp(),
    UNIQUE (database, path)
)
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
    """
    Create the tracking tables if they do not exist yet. A database that predates the
    source_id column needs the one-off ALTER in docs/incremental-parsing.md: it takes
    an exclusive lock on every partition, which is not something to do implicitly in
    the middle of an import.
    """
    with conn.cursor() as cur:
        cur.execute(CREATE_MANIFEST_SQL)
        cur.execute(CREATE_FILES_SQL)
    conn.commit()


def load_file_signatures(conn, database: str) -> ty.Dict[str, str]:
    """The stored ``source path -> signature`` map for a database (empty if none)."""
    ensure_table(conn)
    with conn.cursor() as cur:
        cur.execute(
            f"SELECT path, signature FROM {FILES_TABLE} WHERE database = %s",
            (database,),
        )
        return dict(cur.fetchall())


def resolve_source_ids(
    conn, database: str, paths: ty.Iterable[str]
) -> ty.Dict[str, int]:
    """Ids for these source paths, inserting any the files table has not seen."""
    paths = sorted(set(paths))
    if not paths:
        return {}
    with conn.cursor() as cur:
        execute_values(
            cur,
            f"""
            INSERT INTO {FILES_TABLE} (database, path, signature)
            VALUES %s
            ON CONFLICT (database, path) DO NOTHING
            """,
            [(database, path, "") for path in paths],
            page_size=5000,
        )
        cur.execute(
            f"SELECT path, id FROM {FILES_TABLE} WHERE database = %s AND path = ANY(%s)",
            (database, paths),
        )
        return dict(cur.fetchall())


def store_file_signatures(
    conn,
    database: str,
    signatures: ty.Mapping[str, str],
    dropped: ty.Iterable[str] = (),
) -> None:
    """
    Replace the stored source signatures for a database. Like the record manifest
    this must run only after the load succeeds, so a failed run re-fetches the same
    sources rather than skipping them as unchanged.
    """
    ensure_table(conn)
    dropped = list(dropped)
    with conn.cursor() as cur:
        if dropped:
            # The records go with the source. Leaving them would strand rows whose
            # signature still matches, so a source that came back would be recognised
            # as unchanged and never re-parsed -- while its xrefs stayed retired.
            cur.execute(
                f"""
                DELETE FROM {MANIFEST_TABLE}
                WHERE database = %s
                  AND source_id IN (
                    SELECT id FROM {FILES_TABLE}
                    WHERE database = %s AND path = ANY(%s)
                  )
                """,
                (database, database, dropped),
            )
            cur.execute(
                f"DELETE FROM {FILES_TABLE} WHERE database = %s AND path = ANY(%s)",
                (database, dropped),
            )
        execute_values(
            cur,
            f"""
            INSERT INTO {FILES_TABLE} AS f (database, path, signature)
            VALUES %s
            ON CONFLICT (database, path)
            DO UPDATE SET signature = excluded.signature,
                          updated_at = clock_timestamp()
            WHERE f.signature IS DISTINCT FROM excluded.signature
            """,
            [(database, path, sig) for path, sig in signatures.items()],
            page_size=5000,
        )
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


def load_file_signatures_for(db_url: str, database: str) -> ty.Dict[str, str]:
    """Convenience wrapper that opens its own connection."""
    conn = psycopg2.connect(db_url)
    try:
        return load_file_signatures(conn, database)
    finally:
        conn.close()


def store_signatures(
    conn,
    database: str,
    signatures: ty.Mapping[str, str],
    dropped: ty.Iterable[str] = (),
    sources: ty.Optional[ty.Mapping[str, str]] = None,
) -> None:
    """
    Replace the stored manifest for a database: upsert every current signature and
    remove dropped accessions. Call this only after the database's load succeeds,
    so a failed load leaves the previous manifest intact.

    ``sources`` maps an accession to the source it was read from; it is what lets a
    later import skip an unchanged source without its records looking dropped.
    """
    ensure_partition(conn, database)
    dropped = list(dropped)
    sources = sources or {}
    source_ids = resolve_source_ids(conn, database, sources.values())
    with conn.cursor() as cur:
        if dropped:
            cur.execute(
                f"DELETE FROM {MANIFEST_TABLE} WHERE database = %s AND accession = ANY(%s)",
                (database, dropped),
            )
        # Batched upsert: a per-row INSERT loop is fine for a few thousand HGNC rows
        # but hopeless for ENA's millions, so send them in pages.
        rows = (
            (database, accession, signature, source_ids.get(sources.get(accession)))
            for accession, signature in signatures.items()
        )
        # Without the WHERE, every run rewrites every row -- at ENA scale hundreds of
        # millions of dead tuples for a handful of changed signatures.
        execute_values(
            cur,
            f"""
            INSERT INTO {MANIFEST_TABLE} AS m (database, accession, signature, source_id)
            VALUES %s
            ON CONFLICT (database, accession)
            DO UPDATE SET signature = excluded.signature,
                          source_id = excluded.source_id,
                          updated_at = clock_timestamp()
            WHERE m.signature IS DISTINCT FROM excluded.signature
               OR m.source_id IS DISTINCT FROM excluded.source_id
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
    sources: ty.Optional[ty.Mapping[str, str]] = None,
) -> None:
    """
    Write the two side-channel files a delta parse produces, for the load step:
      * manifest.csv  -- database,accession,signature,source for every current record;
      * deletions.csv -- database,accession for every dropped record.

    The source column is empty for a database that reads one file per import.
    """
    out = Path(output_dir)
    sources = sources or {}
    with (out / MANIFEST_CSV).open("w", newline="") as handle:
        writer = csv.writer(handle)
        for accession, signature in sorted(signatures.items()):
            writer.writerow(
                [database, accession, signature, sources.get(accession, "")]
            )
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
    database's load/release has committed. Reads the (database, accession, signature,
    source) rows from manifest_csv and the dropped (database, accession) rows from
    deletions_csv, then replaces each database's stored signatures.
    """
    signatures: ty.Dict[str, ty.Dict[str, str]] = {}
    sources: ty.Dict[str, ty.Dict[str, str]] = {}
    with Path(manifest_csv).open("r", newline="") as handle:
        for row in csv.reader(handle):
            database, accession, signature = row[:3]
            signatures.setdefault(database, {})[accession] = signature
            source = row[3] if len(row) > 3 else ""
            if source:
                sources.setdefault(database, {})[accession] = source

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
            sources.get(database, {}),
        )


def apply_source_artifacts(conn, sources_csv: ty.Union[str, Path]) -> None:
    """
    Promote the per-source signatures that decide what gets fetched next run. The
    file lists every source seen this run, so a stored source missing from it has
    gone from the input and is forgotten, to be re-fetched if it ever comes back.
    """
    current: ty.Dict[str, ty.Dict[str, str]] = {}
    with Path(sources_csv).open("r", newline="") as handle:
        for database, path, signature in csv.reader(handle):
            current.setdefault(database, {})[path] = signature

    for database, signatures in current.items():
        stored = load_file_signatures(conn, database)
        store_file_signatures(
            conn,
            database,
            signatures,
            dropped=[path for path in stored if path not in signatures],
        )


@attr.s(frozen=True)
class DbDiff:
    """Result of a large-scale diff (see :func:`diff_via_polars`)."""

    to_parse = attr.ib(type=ty.List[str])  # accessions new or changed since last import
    deletions = attr.ib(type=ty.List[str])  # accessions present last import, now absent
    is_bootstrap = attr.ib(type=bool)  # True when there was no prior manifest


def dump_signatures(conn, database: str, path: Path) -> None:
    """
    Stream a database's stored manifest to a headerless ``accession,signature,source``
    CSV.

    This is the only database work an incremental import needs. Keeping it to a
    single COPY of one partition means a diff that lands while the database is busy
    can do nothing but read: no temp tables, no server-side set arithmetic, nothing
    another job can contend with. :func:`diff_via_polars` does the joins.
    """
    ensure_table(conn)
    select = cursor_mogrify(
        conn,
        f"""
        COPY (
            SELECT m.accession, m.signature, coalesce(f.path, '')
            FROM {MANIFEST_TABLE} m
            LEFT JOIN {FILES_TABLE} f ON f.id = m.source_id
            WHERE m.database = %s
        ) TO STDOUT WITH CSV
        """,
        (database,),
    )
    with conn.cursor() as cur, Path(path).open("w") as out:
        cur.copy_expert(select, out)


def cursor_mogrify(conn, statement: str, params: ty.Tuple) -> str:
    """Render a statement with its parameters; COPY cannot take them separately."""
    with conn.cursor() as cur:
        return cur.mogrify(statement, params).decode()


def _scan_signatures(path: Path, columns: ty.List[str]) -> pl.LazyFrame:
    return pl.scan_csv(
        path,
        has_header=False,
        new_columns=columns,
        schema_overrides={column: pl.String for column in columns},
    )


def diff_via_polars(
    stored_csv: Path,
    new_csv: Path,
    scanned_sources: ty.Optional[ty.Collection[str]] = None,
) -> DbDiff:
    """
    Diff this import's signatures against the stored manifest, on the cluster.

    ``stored_csv`` is a headerless ``accession,signature,source`` CSV from
    :func:`dump_signatures` and ``new_csv`` a ``source,accession,signature`` one
    collected from the import; the joins run streaming in polars, so ENA-scale sets never have
    to fit in memory and the database is left alone. An empty ``stored_csv`` is a
    bootstrap: nothing was imported before, so there is nothing to delete and the
    caller parses everything.

    ``scanned_sources`` restricts deletion to records whose source was actually read
    this run. Without it an import that skips an unchanged source would see every one
    of that source's records as absent, and retire the lot.

    ENA can emit the same location-based accession twice in one snapshot, so the new
    side is deduplicated before the join.
    """
    if Path(stored_csv).stat().st_size == 0:
        return DbDiff(to_parse=[], deletions=[], is_bootstrap=True)

    if Path(new_csv).stat().st_size == 0:
        # Every source was skipped as unchanged. Nothing to parse, but a source that
        # has gone still has records to retire, so the anti-join must still run.
        new = pl.LazyFrame(schema={"accession": pl.String, "signature": pl.String})
    else:
        new = (
            _scan_signatures(new_csv, ["source", "accession", "signature"])
            .select("accession", "signature")
            .unique(subset="accession", keep="any")
        )
    stored = _scan_signatures(stored_csv, ["accession", "signature", "source"])

    to_parse = (
        new.join(stored, on="accession", how="left", suffix="_stored")
        .filter(
            pl.col("signature_stored").is_null()
            | (pl.col("signature_stored") != pl.col("signature"))
        )
        .select("accession")
    )
    deletions = stored
    if scanned_sources is not None:
        deletions = deletions.filter(pl.col("source").is_in(list(scanned_sources)))
    deletions = deletions.join(new, on="accession", how="anti").select("accession")

    return DbDiff(
        to_parse=to_parse.collect(engine="streaming")["accession"].to_list(),
        deletions=deletions.collect(engine="streaming")["accession"].to_list(),
        is_bootstrap=False,
    )
