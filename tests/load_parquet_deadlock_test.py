# -*- coding: utf-8 -*-

"""
``bin/load-parquet`` talks to Postgres over two connections: DuckDB's
``postgres`` attachment does the TRUNCATE and the INSERT, and ``_run_post_load``
opens its own psycopg2 connection for the merge.

The attachment keeps its Postgres transaction open until the DuckDB connection
closes, so it is still holding the ACCESS EXCLUSIVE lock TRUNCATE took. Running
the post-load before closing it deadlocks the merge's CREATE INDEX against a
lock the same process holds — no timeout, no error, just a session blocked on
itself. The rfam load sat that way for 18 hours before anyone noticed.
"""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOAD_PARQUET = ROOT / "bin" / "load-parquet"


def _main_body():
    tree = ast.parse(LOAD_PARQUET.read_text())
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "main":
            return node
    raise AssertionError("bin/load-parquet has no main()")


def _call_lines(fn, predicate):
    return [
        node.lineno
        for node in ast.walk(fn)
        if isinstance(node, ast.Call) and predicate(node.func)
    ]


def test_duckdb_connection_closed_before_post_load():
    fn = _main_body()

    closes = _call_lines(
        fn,
        lambda f: isinstance(f, ast.Attribute)
        and f.attr == "close"
        and isinstance(f.value, ast.Name)
        and f.value.id == "con",
    )
    post_loads = _call_lines(
        fn, lambda f: isinstance(f, ast.Name) and f.id == "_run_post_load"
    )

    assert closes, "main() never closes the DuckDB connection"
    assert post_loads, "main() never calls _run_post_load"
    assert min(closes) < min(post_loads), (
        "con.close() must come before _run_post_load: the DuckDB attachment "
        "holds TRUNCATE's ACCESS EXCLUSIVE lock until it closes, so the "
        "post-load's own connection deadlocks against it"
    )
