import asyncio
import datetime
import importlib.util
from pathlib import Path

import polars as pl

BIN = Path(__file__).resolve().parents[2] / "bin" / "litscan-search-job.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("litscan_search_job", BIN)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


search_job = _load_script()


def test_hit_after_many_zero_hit_searches(monkeypatch):
    # polars infers column types from the first 100 rows; all-null pmcids
    # there made the first real hit panic on concatenation.
    hit = "urs0000000001"

    async def fake_search(session, limiter, semaphore, job_id, date, search_limit):
        if job_id == hit:
            return {
                "hit_count": [1],
                "pmcids": ["PMC123"],
                "cite_counts": [4],
                "status": ["success"],
            }
        return {
            "hit_count": [0],
            "pmcids": [None],
            "cite_counts": [None],
            "status": ["success"],
        }

    monkeypatch.setattr(search_job, "search_article_async", fake_search)
    job_ids = [f"none_{i}" for i in range(100)] + [hit]
    df = pl.DataFrame(
        {"job_id": job_ids, "finished": [datetime.datetime(2026, 1, 1)] * 101}
    )

    result = asyncio.run(search_job.fetch_all_epmc_data(df, search_limit=10))

    row = result.filter(pl.col("job_id") == hit)
    assert row["pmcids"].to_list() == [["PMC123"]]
    assert row["cite_counts"].to_list() == [[4]]
