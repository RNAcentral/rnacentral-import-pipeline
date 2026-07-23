# -*- coding: utf-8 -*-

"""
Shared machinery for incremental ("delta") parsing.

A database's import can skip the expensive per-record parsing for records that
have not changed since the last import. To decide what changed we keep a manifest
of ``accession -> signature`` (a hash of the raw record) from the last successful
import, diff the new input against it, and only fully parse the new/changed
records. Records that vanished from the input become an explicit deletion list.

This module is database-agnostic: signature/diff are pure, and the manifest table
is a single generic ``rnc_import_manifest`` keyed by (database, accession). See
docs/incremental-parsing.md.
"""

import csv
import hashlib
import json
from pathlib import Path
import typing as ty

import attr
import psycopg2

MANIFEST_TABLE = "rnacen.rnc_import_manifest"

MANIFEST_CSV = "manifest.csv"
DELETIONS_CSV = "deletions.csv"

CREATE_MANIFEST_SQL = f"""
CREATE TABLE IF NOT EXISTS {MANIFEST_TABLE} (
    database   text        NOT NULL,
    accession  text        NOT NULL,
    signature  text        NOT NULL,
    updated_at timestamptz NOT NULL DEFAULT clock_timestamp(),
    PRIMARY KEY (database, accession)
)
"""


def record_signature(raw: ty.Any) -> str:
    """
    A stable hash of a raw record. Canonical JSON (sorted keys) so the signature
    depends only on content, not key order or whitespace. The whole record is
    hashed, so any change to any field is caught.
    """
    canonical = json.dumps(raw, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


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
    """Create the manifest table if it does not exist yet."""
    with conn.cursor() as cur:
        cur.execute(CREATE_MANIFEST_SQL)
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
    ensure_table(conn)
    dropped = list(dropped)
    with conn.cursor() as cur:
        if dropped:
            cur.execute(
                f"DELETE FROM {MANIFEST_TABLE} WHERE database = %s AND accession = ANY(%s)",
                (database, dropped),
            )
        for accession, signature in signatures.items():
            cur.execute(
                f"""
                INSERT INTO {MANIFEST_TABLE} (database, accession, signature)
                VALUES (%s, %s, %s)
                ON CONFLICT (database, accession)
                DO UPDATE SET signature = excluded.signature,
                              updated_at = clock_timestamp()
                """,
                (database, accession, signature),
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
    Promote a delta parse's manifest into rnc_import_manifest. Call only after the
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
