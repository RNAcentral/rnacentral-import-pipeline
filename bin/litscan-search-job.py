#!/usr/bin/env python3
import argparse
import asyncio
import datetime
import logging
import os
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import ParseError

import aiohttp
import polars as pl
import psycopg2
import requests
from aiolimiter import AsyncLimiter

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"


def search_article(
    job_id: str, date: datetime.date, search_limit: int
) -> dict[str, list]:
    page = "*"
    articles_list = []
    search_date = (
        f" AND FIRST_PDATE:[{date.strftime('%Y-%m-%d')} TO {datetime.date.today().strftime('%Y-%m-%d')}]"
        if date
        else ""
    )

    while len(articles_list) < search_limit and page:
        query = (
            f'search?query=("{job_id}" AND IN_EPMC:Y AND OPEN_ACCESS:Y AND NOT SRC:PPR{search_date})'
            f"&pageSize=500&cursorMark={page}"
        )
        print(query)
        try:
            response = requests.get(EUROPE_PMC + query, timeout=60)
            response.raise_for_status()
            articles = response.text
        except requests.exceptions.RequestException as e:
            logger.warning("Error fetching article list from Europe PMC: %s", e)
            return {"pmcids": [], "cite_counts": []}
        print(articles)
        try:
            root = ET.fromstring(articles)
        except ParseError:
            logger.warning(f"parse failure for {job_id}")
            logger.warning(articles)
            return {"pmcids": [], "cite_counts": []}

        articles_list.extend(
            [
                (item.find("pmcid").text, int(item.find("citedByCount").text))
                for item in root.findall("./resultList/result")
                if item.find("pmcid") is not None
                and item.find("citedByCount") is not None
            ]
        )

        try:
            next_page = root.find("nextCursorMark").text
        except AttributeError:
            next_page = None

        page = next_page

    if not articles_list:
        logger.warning(f"No articles found for {job_id}")
        return {"pmcids": [], "cite_counts": []}

    pmcids, hit_counts = zip(*articles_list)

    return {"pmcids": list(pmcids), "cite_counts": list(hit_counts)}


async def search_article_async(
    session: aiohttp.ClientSession,
    limiter: AsyncLimiter,
    job_id: str,
    date: datetime.date,
    search_limit: int,
) -> dict[str, list]:
    page = "*"
    articles_list = []
    search_date = (
        f" AND FIRST_PDATE:[{date.strftime('%Y-%m-%d')} TO {datetime.date.today().strftime('%Y-%m-%d')}]"
        if date
        else ""
    )

    while len(articles_list) < search_limit and page:
        query = (
            f'search?query=("{job_id}" AND IN_EPMC:Y AND OPEN_ACCESS:Y AND NOT SRC:PPR{search_date})'
            f"&pageSize=500&cursorMark={page}"
        )

        max_retries = 3
        articles = None
        for attempt in range(max_retries):
            try:
                async with limiter:
                    async with session.get(EUROPE_PMC + query, timeout=60) as response:
                        if (
                            response.status in (429, 500, 502, 503, 504)
                            and attempt < max_retries - 1
                        ):
                            delay = 2**attempt
                            logger.warning(
                                "Got %d for %s, retrying in %ds (attempt %d/%d)",
                                response.status,
                                job_id,
                                delay,
                                attempt + 1,
                                max_retries,
                            )
                            await asyncio.sleep(delay)
                            continue
                        response.raise_for_status()
                        articles = await response.text()
            except aiohttp.ClientError as e:
                if attempt < max_retries - 1:
                    delay = 2**attempt
                    logger.warning(
                        "Network error for %s: %s, retrying in %ds", job_id, e, delay
                    )
                    await asyncio.sleep(delay)
                    continue
                logger.warning(
                    "Network error for %s after %d attempts: %s", job_id, max_retries, e
                )
                return {"hit_count": [], "pmcids": [], "cite_counts": []}
            except asyncio.TimeoutError:
                if attempt < max_retries - 1:
                    delay = 2**attempt
                    logger.warning("Timeout for %s, retrying in %ds", job_id, delay)
                    await asyncio.sleep(delay)
                    continue
                logger.warning("Timeout for %s after %d attempts", job_id, max_retries)
                return {"hit_count": [], "pmcids": [], "cite_counts": []}
            else:
                break
        if articles is None:
            return {"hit_count": [], "pmcids": [], "cite_counts": []}

        try:
            root = ET.fromstring(articles)
        except ParseError:
            logger.warning(f"Parse failure for {job_id}")
            return {"hit_count": [], "pmcids": [], "cite_counts": []}
        hit_count_el = root.find("hitCount")
        if hit_count_el is None:
            logger.warning("No hitCount in response for %s, skipping page", job_id)
            break
        hit_count = int(hit_count_el.text)

        articles_list.extend(
            [
                (
                    hit_count,
                    item.find("pmcid").text,
                    int(item.find("citedByCount").text),
                )
                for item in root.findall("./resultList/result")
                if item.find("pmcid") is not None
                and item.find("citedByCount") is not None
            ]
        )

        try:
            next_page = root.find("nextCursorMark").text
        except AttributeError:
            next_page = None

        page = next_page

    if not articles_list:
        return {"hit_count": [], "pmcids": [], "cite_counts": []}

    hit_counts, pmcids, cite_counts = zip(*articles_list)

    return {
        "hit_count": list(hit_counts),
        "pmcids": list(pmcids),
        "cite_counts": list(cite_counts),
    }


async def fetch_all_epmc_data(
    df: pl.DataFrame, search_limit: int = 100
) -> pl.DataFrame:
    """
    Takes a Polars DataFrame with 'job_id' and 'date' columns,
    runs the async Europe PMC search, and returns a new DataFrame with the results attached.
    """
    # 10 requests per 1 second
    limiter = AsyncLimiter(10, 1.0)

    # Extract data to Python lists
    job_ids = df["job_id"].to_list()
    dates = df["finished"].to_list()

    # Open a single persistent HTTP session for all requests
    async with aiohttp.ClientSession() as session:
        # Create a list of tasks (they don't run yet)
        tasks = [
            search_article_async(session, limiter, jid, d, search_limit)
            for jid, d in zip(job_ids, dates)
        ]

        # Run them all concurrently and wait for them to finish
        print(f"Executing {len(tasks)} requests at 10 req/sec...")
        results = await asyncio.gather(*tasks)

    results_df = pl.DataFrame(results)
    final_df = pl.concat([df, results_df], how="horizontal")

    return final_df


def articles_list(date: pl.DataFrame, search_limit: int):
    search_results = asyncio.run(fetch_all_epmc_data(date, search_limit=search_limit))

    return search_results


def get_search_dates(conn, job_ids, chunk_size=5000):
    """
    PostgreSQL/psycopg2 optimized version using the ANY() operator.
    """
    query = "SELECT job_id, finished FROM litscan_job WHERE job_id = ANY(%s)"

    df = pl.read_database(
        query=query,
        connection=conn,
        execute_options={"parameters": (job_ids,)},
        schema_overrides={"finished": pl.Datetime("us")},
    )
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Search many RNA IDs in Europe PMC and store results in CSV"
    )
    parser.add_argument(
        "--job-id", required=True, help="Text file with one job_id per line"
    )
    parser.add_argument(
        "--output", required=True, help="Filename to write search results to"
    )
    args = parser.parse_args()

    job_ids = [jid.lower() for jid in Path(args.job_id).read_text().split()]

    conn_str = os.environ["PSYCOPG_CONN"]
    conn = psycopg2.connect(conn_str)
    with conn.cursor() as cur:
        last_search = get_search_dates(cur, job_ids)
    conn.close()

    search_results = articles_list(last_search, search_limit=1_000_000)

    search_results = search_results.explode("pmcids", "cite_counts", "hit_count")

    search_results.sort(by="pmcids").write_csv(args.output, include_header=False)


if __name__ == "__main__":
    main()
