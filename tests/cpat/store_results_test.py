# -*- coding: utf-8 -*-

"""
Sequences with no ORF are written with null scores, but the live load_cpat was
created before that with NOT NULL on every score column and create_load.sql is
never re-applied. The pgloader ctl drops the constraints in BEFORE LOAD DO; the
parquet branch of store_results must do the same before load-parquet runs, or
the first no-ORF row kills the load.
"""

import re
from pathlib import Path

from rnacentral_pipeline import schemas

WORKFLOW = Path(__file__).resolve().parents[2] / "workflows" / "cpat.nf"


def test_parquet_branch_drops_score_not_nulls_before_loading():
    text = WORKFLOW.read_text()
    load = text.index("load-parquet load_cpat ")
    alter = re.search(r"ALTER TABLE load_cpat[^\n]*", text[:load])
    assert alter, "no ALTER TABLE load_cpat before load-parquet"
    nullable = [f.name for f in schemas.CPAT_RESULTS if f.nullable]
    for column in nullable:
        assert f"ALTER COLUMN {column} DROP NOT NULL" in alter.group(0), column
