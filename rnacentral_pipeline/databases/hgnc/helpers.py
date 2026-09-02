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

import collections as coll
import hashlib
import json
import logging
import threading
import typing as ty
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from rnacentral_pipeline.databases.hgnc.data import HgncEntry

LOGGER = logging.getLogger(__name__)

ENSEMBL_SEQUENCE_URL = "https://rest.ensembl.org/sequence/id"

# Ensembl's POST endpoint accepts at most 50 ids per request.
ENSEMBL_BATCH_SIZE = 50
ENSEMBL_TIMEOUT = 120
# Each request takes ~12s and covers 50 genes, so this is well under Ensembl's
# 15 req/s allowance. Fetching is entirely network bound, so threads are enough.
ENSEMBL_WORKERS = 4

REFSEQ_DBID = 9
GTRNADB_DBID = 8
HUMAN_TAXID = 9606


def url(entry: HgncEntry) -> str:
    return (
        f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{entry.hgnc_id}"
    )


def load(path: Path) -> ty.List[HgncEntry]:
    data = []
    with path.open("r") as handle:
        raw_data = json.load(handle)
        for raw in raw_data["response"]["docs"]:
            data.append(HgncEntry.from_raw(raw))
    return data


def description(entry: HgncEntry) -> str:
    return f"Homo sapiens (human) {entry.name}"


def gene(entry: HgncEntry) -> str:
    return entry.symbol


def gene_synonyms(hgnc: HgncEntry) -> ty.List[str]:
    return hgnc.previous_symbols


def md5(sequence: str) -> str:
    sequence = sequence.replace("U", "T").upper()
    m = hashlib.md5(sequence.encode())
    return m.hexdigest()


def _longest(rows) -> ty.Dict[str, str]:
    """
    Collapse (key, urs, length) rows to {key: urs}, keeping the longest
    sequence for each key.
    """
    best: ty.Dict[str, ty.Tuple[str, int]] = {}
    for key, urs, length in rows:
        if key not in best or length > best[key][1]:
            best[key] = (urs, length)
    return {key: urs for key, (urs, _) in best.items()}


def refseq_mapping(conn, refseq_ids: ty.List[str]) -> ty.Dict[str, str]:
    """
    Map every RefSeq id to the longest URS it corresponds to.

    Scoped to the RefSeq xref partition. Without the dbid the planner probes
    all ~180 xref partitions, which turns this from seconds into minutes.
    RefSeq is the only database that resolves any of these ids anyway.
    """
    if not refseq_ids:
        return {}

    ids = sorted(set(refseq_ids))
    query = """
    select k, rna.urs, rna.len
    from (
        select acc.parent_ac as k, acc.accession
        from rnc_accessions acc where acc.parent_ac = ANY(%s)
      union all
        select acc.external_id, acc.accession
        from rnc_accessions acc where acc.external_id = ANY(%s)
      union all
        select acc.optional_id, acc.accession
        from rnc_accessions acc where acc.optional_id = ANY(%s)
    ) m
    join xref on xref.ac = m.accession
    join rna on rna.urs = xref.urs
    where xref.taxid = %s and xref.deleted = 'N' and xref.dbid = %s
    """
    with conn.cursor() as cur:
        cur.execute(query, (ids, ids, ids, HUMAN_TAXID, REFSEQ_DBID))
        return _longest(cur)


def gtrnadb_mapping(conn, gtrnadb_ids: ty.List[str]) -> ty.Dict[str, str]:
    """
    Map every gtRNAdb id to the longest URS it corresponds to.
    """
    if not gtrnadb_ids:
        return {}

    query = """
    select acc.optional_id, rna.urs, rna.len
    from xref
    join rnc_accessions acc on acc.accession = xref.ac
    join rna on rna.urs = xref.urs
    where xref.taxid = %s and xref.deleted = 'N' and xref.dbid = %s
      and acc.database = 'GTRNADB' and acc.optional_id = ANY(%s)
    """
    with conn.cursor() as cur:
        cur.execute(query, (HUMAN_TAXID, GTRNADB_DBID, sorted(set(gtrnadb_ids))))
        return _longest(cur)


_LOCAL = threading.local()


def _ensembl_session() -> requests.Session:
    """
    One session per thread, as requests sessions are not thread safe.

    Ensembl REST returns intermittent 5xx and rate limits with 429. A batch is
    50 genes, so an unretried failure loses all of them silently.
    """
    session = getattr(_LOCAL, "session", None)
    if session is None:
        session = requests.Session()
        retry = Retry(
            total=5,
            backoff_factor=2,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["POST"],
            respect_retry_after_header=True,
        )
        session.mount("https://", HTTPAdapter(max_retries=retry))
        _LOCAL.session = session
    return session


def _fetch_batch(batch: ty.List[str]) -> ty.List[dict]:
    payload = {"ids": batch, "type": "cdna", "multiple_sequences": 1}
    try:
        response = _ensembl_session().post(
            ENSEMBL_SEQUENCE_URL,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            json=payload,
            timeout=ENSEMBL_TIMEOUT,
        )
        response.raise_for_status()
    except requests.exceptions.RequestException as err:
        LOGGER.error("Gave up on a batch, losing %i genes: %s", len(batch), err)
        return []
    return response.json()


def ensembl_cdna(gene_ids: ty.List[str]) -> ty.Dict[str, ty.List[str]]:
    """
    Fetch the spliced transcript sequences for each Ensembl gene, batched.

    We ask for cdna rather than the default genomic sequence because RNAcentral
    stores spliced transcripts; the genomic sequence of a multi-exon gene can
    never md5 match one of ours.
    """
    found: ty.Dict[str, ty.List[str]] = coll.defaultdict(list)
    if not gene_ids:
        return found

    ids = sorted(set(gene_ids))
    batches = [
        ids[index : index + ENSEMBL_BATCH_SIZE]
        for index in range(0, len(ids), ENSEMBL_BATCH_SIZE)
    ]
    done = 0
    with ThreadPoolExecutor(max_workers=ENSEMBL_WORKERS) as pool:
        for records in pool.map(_fetch_batch, batches):
            # Ensembl silently drops ids it cannot resolve, so match on the
            # echoed query rather than assuming the response lines up with the
            # request.
            for record in records:
                found[record.get("query", record["id"])].append(record["seq"])
            done += 1
            LOGGER.info("Fetched Ensembl batch %i/%i", done, len(batches))
    return found


def ensembl_mapping(conn, gene_ids: ty.List[str]) -> ty.Dict[str, str]:
    """
    Map Ensembl gene ids to URS by md5 matching each of their transcripts.
    """
    sequences = ensembl_cdna(gene_ids)
    if not sequences:
        return {}

    digests: ty.Dict[str, ty.List[str]] = coll.defaultdict(list)
    for gene_id, seqs in sequences.items():
        for seq in seqs:
            digests[md5(seq)].append(gene_id)

    with conn.cursor() as cur:
        cur.execute(
            "select md5, urs, len from rna where md5 = ANY(%s)", (list(digests),)
        )
        rows = [
            (gene_id, upi, length)
            for digest, upi, length in cur
            for gene_id in digests[digest]
        ]
    return _longest(rows)


def sequence_mapping(conn, urs_ids: ty.List[str]) -> ty.Dict[str, str]:
    """
    Fetch the sequence for every URS we resolved.
    """
    if not urs_ids:
        return {}

    query = """
    select urs, coalesce(seq_short, seq_long)
    from rna where urs = ANY(%s)
    """
    with conn.cursor() as cur:
        cur.execute(query, (sorted(set(urs_ids)),))
        return dict(cur)


SO_TERMS = {
    "RNA, long non-coding": "SO:0001877",
    "RNA, Y": "SO:0000405",
    "RNA, cluster": "SO:0000655",
    "RNA, micro": "SO:0000276",
    "RNA, misc": "SO:0000655",
    "RNA, small nuclear": "SO:0000274",
    "RNA, small nucleolar": "SO:0000275",
    "RNA, transfer": "SO:0000253",
    "RNA, vault": "SO:0000404",
}


def so_term(entry: HgncEntry) -> ty.Optional[str]:
    """
    Get the SO term for an entry, or None if we do not know the RNA type.

    Returning None rather than raising means a newly introduced HGNC locus
    type costs us one entry instead of the whole run.
    """
    if entry.hgnc_rna_type == "RNA, ribosomal":
        if "5S ribosomal RNAs" in entry.gene_groups:
            return "SO:0000652"
        if "12S RNA" in entry.gene_groups:
            return "SO:0002344"
        return "SO:0000252"

    term = SO_TERMS.get(entry.hgnc_rna_type)
    if term is None:
        LOGGER.warning(
            "Unknown type of RNA (%s) for %s", entry.hgnc_rna_type, entry.hgnc_id
        )
    return term
