# -*- coding: utf-8 -*-

"""
Delta loading for the unified Ensembl import: a geneset whose files have not changed
is neither fetched nor parsed, and its xrefs are left alone. See
docs/incremental-parsing.md; the source machinery is the one ENA uses.
"""

import contextlib
import csv
import hashlib
import logging
import time
import typing as ty
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from rnacentral_pipeline.databases import manifest

LOGGER = logging.getLogger(__name__)

# One import yields all of these, and each needs its own manifest rows: a database with
# none is loaded FULL, which retires everything a partial parse did not re-emit.
DATABASES = ("ENSEMBL", "ENSEMBL_GENCODE", "TAIR")

Row = ty.Tuple[str, ...]


def mirror(ftp_url: str, nfs_path: str) -> ty.Callable[[str], str]:
    """Reads an FTP url from its copy on /nfs; the FTP takes ~19 s per file from a compute node."""
    return lambda url: Path(url.replace(ftp_url, nfs_path, 1)).read_text()


def source_signature(
    embl_url: str, gff_url: str, fetch: ty.Callable[[str], str]
) -> str:
    """
    A hash of the two files we read, from the md5sum.txt in their directory. If that
    cannot be read the signature carries the time, so it never matches a stored one
    and the geneset is parsed again: slower, never wrong.
    """
    names = {u.rsplit("/", 1)[1]: u for u in (embl_url, gff_url)}
    try:
        listing = fetch(embl_url.rsplit("/", 1)[0] + "/md5sum.txt")
        md5s = {}
        for line in listing.splitlines():
            checksum, _, name = line.partition(" ")
            md5s[name.strip().lstrip("*").rsplit("/", 1)[-1]] = checksum
        parts = [f"{url}\t{md5s[name]}" for name, url in sorted(names.items())]
    except (OSError, KeyError) as err:
        LOGGER.warning("No checksums for %s (%s); it will be parsed", embl_url, err)
        return f"unverified-{time.time()}"
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()


def signatures(
    rows: ty.Iterable[Row], read: ty.Callable[[str], str]
) -> ty.Dict[str, str]:
    """Source signature for each geneset row, keyed by its EMBL url."""
    rows = list(rows)
    with ThreadPoolExecutor(max_workers=16) as pool:
        sigs = pool.map(lambda r: source_signature(r[2], r[3], read), rows)
        return {row[2]: sig for row, sig in zip(rows, sigs)}


def select(rows: ty.Iterable[Row], species: ty.Collection[str]) -> ty.List[Row]:
    """The rows for the named species, or all of them when none are named."""
    chosen = {s.lower() for s in species}
    return [r for r in rows if not chosen or r[0].lower() in chosen]


def plan(
    rows: ty.Sequence[Row],
    picked: ty.Collection[Row],
    stored: ty.Mapping[str, str],
    current: ty.Mapping[str, str],
) -> ty.Tuple[ty.List[Row], ty.List[str], ty.Dict[str, str]]:
    """
    Split the genesets into those to fetch and the ones that have gone.

    ``picked`` are the rows to look at this run and ``current`` their signatures.
    A geneset that was not picked keeps the signature it had, so a partial run cannot
    mark it up to date, and is not counted as gone.
    """
    picked_urls = {r[2] for r in picked}
    to_fetch = [
        r + (current[r[2]],) for r in picked if stored.get(r[2]) != current[r[2]]
    ]
    sources = {r[2]: current[r[2]] for r in picked}
    sources.update(
        {r[2]: stored[r[2]] for r in rows if r[2] not in picked_urls and r[2] in stored}
    )
    gone = sorted(set(stored) - {r[2] for r in rows})
    return to_fetch, gone, sources


def deletions(
    conn, manifest_csv: Path, scanned: ty.Collection[str], workdir: Path
) -> ty.List[Row]:
    """
    The (database, accession) rows stored for a geneset that was read this run but
    are not in what it now yields.
    """
    with contextlib.ExitStack() as stack, manifest_csv.open(newline="") as handle:
        writers = {
            db: csv.writer(
                stack.enter_context((workdir / f"new-{db}.csv").open("w", newline=""))
            )
            for db in DATABASES
        }
        for database, accession, signature, source in csv.reader(handle):
            writers[database].writerow((source, accession, signature))

    found = []
    for database in DATABASES:
        stored = workdir / f"stored-{database}.csv"
        manifest.dump_signatures(conn, database, stored)
        current = workdir / f"new-{database}.csv"
        result = manifest.diff_manifests(stored, current, scanned)
        found.extend((database, accession) for accession in result.deletions)
    return found
