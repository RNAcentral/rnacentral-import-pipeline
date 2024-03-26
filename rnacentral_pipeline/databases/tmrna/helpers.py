import asyncio
import io
import typing as ty

import aiohttp
from Bio import Entrez, SeqIO
from more_itertools import chunked
from throttler import throttle

Entrez.email = "rnacentral@gmail.com"


@throttle(rate_limit=2, period=1.0)
async def fetch_records(session, accessions: ty.List[str]):
    try:
        accession_str = ",".join(accessions)
        async with session.get(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession_str}&rettype=gb&retmode=text"
        ) as response:
            records_text = await response.text()
            handle = io.StringIO(records_text)
            taxids = {}
            for record in SeqIO.parse(handle, format="genbank"):
                taxids[record.id] = extract_taxid(record)
            return taxids
    except Exception as e:
        print(f"Error fetching records: {e}")
        return {accession: None for accession in accessions}


def extract_taxid(record):
    for feature in record.features:
        if feature.type == "source" and "db_xref" in feature.qualifiers:
            db_xrefs = feature.qualifiers["db_xref"]
            for db_xref in db_xrefs:
                if "taxon" in db_xref:
                    return int(db_xref.split(":")[1])
    return None


async def get_taxids_for_accessions(accessions, batch_size=100, requests_per_second=3):
    async with aiohttp.ClientSession() as session:
        tasks = []
        for batch in chunked(accessions, batch_size):
            tasks.append(fetch_records(session, batch))

        taxids = {}
        for result in asyncio.as_completed(tasks):
            mapping = await result
            for accession, taxid in mapping.items():
                if taxid is None:
                    continue
                taxids[accession] = taxid
    return taxids


def taxid_mapping(
    rows: ty.List[ty.Dict[str, str]], batch_size=100
) -> ty.Dict[str, int]:
    accessions = []
    for row in rows:
        acc = row["Instances"].split(",")
        accessions.extend(a.split("/")[0] for a in acc)
    return asyncio.run(get_taxids_for_accessions(accessions, batch_size=batch_size))
