# -*- coding: utf-8 -*-

"""
Delta (incremental) parsing support for ENA.

ENA is fetched as EMBL flat-file ``.ncr`` chunks and normally every record in every
chunk is run through ribotyper and the full parse, even when almost nothing changed
between snapshots. This module provides the two cheap passes that let the workflow
skip the expensive work for unchanged records:

* :func:`iter_signatures` — ``accession -> sha256(raw record)`` for every record in a
  chunk, with no ribotyper and no database access;
* :func:`filter_records` — copy through only the records whose accession is in a
  supplied set, so ribotyper and the parse run on the changed subset alone.

Both derive the accession with the *same* ``SeqIO`` ``record.id`` the real parser
uses (see :mod:`rnacentral_pipeline.databases.ena.parser` and
``ena/helpers.accession``), so the signature keys, the to-parse filter and the
explicit deletion list all agree byte-for-byte. See docs/incremental-parsing-ena.md.
"""

import io
import typing as ty
from pathlib import Path

from Bio import SeqIO

from rnacentral_pipeline.databases import manifest

# The sentinel a bootstrap diff (no prior manifest) writes as the sole line of the
# to-parse file, meaning "keep every record" -- avoids materialising a to-parse set
# of every accession on the first delta run. Kept in sync with ena_delta_diff.
KEEP_ALL = "__ALL__"


def _record_id(block: str) -> ty.Optional[str]:
    """
    The parser's ``record.id`` for a single raw EMBL record block, or None if the
    block does not parse as a record (blank trailing text). Uses SeqIO so the id is
    identical to what parser.py sees -- never a hand-rolled accession scan.
    """
    for record in SeqIO.parse(io.StringIO(block), "embl"):
        return record.id
    return None


def iter_raw_records(path: Path) -> ty.Iterator[ty.Tuple[str, str]]:
    """
    Yield ``(record_id, raw_block)`` for each EMBL record in a ``.ncr`` file.

    Records are delimited exactly as ``utils/split-ena`` delimits them: a record ends
    at a line beginning ``//``. Blocks that do not parse to a record (e.g. trailing
    whitespace after the final ``//``) are skipped.
    """
    block: ty.List[str] = []
    with Path(path).open("r") as handle:
        for line in handle:
            block.append(line)
            if line.startswith("//"):
                raw = "".join(block)
                block = []
                record_id = _record_id(raw)
                if record_id is not None:
                    yield record_id, raw
    if block:
        raw = "".join(block)
        record_id = _record_id(raw)
        if record_id is not None:
            yield record_id, raw


def iter_signatures(path: Path) -> ty.Iterator[ty.Tuple[str, str]]:
    """Yield ``(accession, signature)`` for every record in an ENA chunk."""
    for record_id, raw in iter_raw_records(path):
        yield record_id, manifest.text_signature(raw)


def write_signatures(path: Path, output: ty.IO[str]) -> int:
    """
    Write ``accession,signature`` CSV rows for a chunk to ``output``. Returns the
    number of records written. No database column here -- the diff step adds it once,
    so the per-chunk files stay small.
    """
    import csv

    writer = csv.writer(output)
    count = 0
    for accession, signature in iter_signatures(path):
        writer.writerow([accession, signature])
        count += 1
    return count


def _load_keep_set(only: Path) -> ty.Optional[ty.Set[str]]:
    """
    Read a to-parse file into a set of accessions, or None when it is the KEEP_ALL
    sentinel (bootstrap: keep everything).
    """
    accessions: ty.Set[str] = set()
    with Path(only).open("r") as handle:
        for line in handle:
            accession = line.strip()
            if not accession:
                continue
            if accession == KEEP_ALL:
                return None
            accessions.add(accession)
    return accessions


def filter_records(path: Path, only: Path, output: Path) -> int:
    """
    Copy records from ``path`` to ``output`` keeping only those whose accession is
    listed in ``only`` (the to-parse file). The KEEP_ALL sentinel copies everything.

    Returns the number of records written. Zero means the whole chunk was unchanged;
    the workflow then skips ribotyper and the parse for that chunk entirely.
    """
    keep = _load_keep_set(only)
    written = 0
    with Path(output).open("w") as out:
        for record_id, raw in iter_raw_records(path):
            if keep is None or record_id in keep:
                out.write(raw)
                written += 1
    return written
