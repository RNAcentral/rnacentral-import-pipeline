from rnacentral_pipeline.rnacentral.release import run


class FakeCursor:
    def __init__(self):
        self.calls = []
        self._last = ""

    def execute(self, sql, params=None):
        self.calls.append((sql, params))
        self._last = sql

    def executemany(self, sql, seq_of_params):
        # _record_pending_stats batch-inserts the loaded dbids this way.
        for params in seq_of_params:
            self.execute(sql, params)

    def fetchall(self):
        # functions.apply asks which functions were already applied -- pretend none,
        # so it proceeds to (re)deploy them. It expects (schema, name, sha) triples,
        # so handing it the release tuples would blow up in _applied_shas.
        if "applied_functions" in self._last:
            return []
        # run() asks which releases are pending loading (status = 'L').
        if "FROM rnacen.rnc_release" in self._last:
            return [(9, 123)]
        return []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


class FakeConnection:
    def __init__(self):
        self.cursor_obj = FakeCursor()
        self.commits = 0
        self.autocommit = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def cursor(self):
        return self.cursor_obj

    def commit(self):
        self.commits += 1

    def rollback(self):
        pass


def _run_capture(monkeypatch, **kwargs):
    conn = FakeConnection()
    # _connect passes keepalive/options kwargs, so the stub must accept them.
    monkeypatch.setattr(run.psycopg2, "connect", lambda *a, **k: conn)
    run.run("postgres://example", **kwargs)
    return conn.cursor_obj.calls


def test_run_patches_functions_and_checks_once(monkeypatch):
    sql_calls = [sql for sql, _ in _run_capture(monkeypatch)]

    # Functions are deployed from database_functions/ in (schema, name) order,
    # interleaved with tracking-table INSERTs, so match on content rather than
    # position.
    # do_pel_exchange is patched for declarative partitioning support.
    assert any(
        "CREATE OR REPLACE FUNCTION rnc_load_xref.do_pel_exchange" in s
        for s in sql_calls
    )

    # load_xref is patched to drop the per-database do_checks call.
    load_xref = next(
        s
        for s in sql_calls
        if "CREATE OR REPLACE FUNCTION rnc_load_xref.load_xref" in s
    )
    assert "perform rnc_load_xref.do_checks" not in load_xref

    # The xref-join functions are patched to add a partition-pruning predicate.
    pp1 = next(
        s for s in sql_calls if "FUNCTION rnc_load_xref.populate_pel_tables1" in s
    )
    pp2 = next(
        s for s in sql_calls if "FUNCTION rnc_load_xref.populate_pel_tables2" in s
    )
    lumv = next(
        s
        for s in sql_calls
        if "FUNCTION rnc_load_xref.load_upi_max_versions_table" in s
    )
    assert "x.dbid = v_dbid" in pp1
    assert "x.dbid = v_dbid" in pp2
    assert "PREVIOUS_XREF.DBID = p_in_dbid" in lumv

    # log_release_end_atx is patched to run its stats query through EXECUTE
    # format(), so the dbid arrives as a literal and the xref partitions prune.
    # The old correlated EXISTS subqueries (p.dbid = x.dbid) are gone entirely.
    lrea = next(s for s in sql_calls if "FUNCTION rnc_logging.log_release_end_atx" in s)
    assert "x.dbid = %1$s" in lrea
    assert "p.dbid    = x.dbid" not in lrea

    # The load_md5_new_sequences index is created to support set_comparable_prot_upi.
    assert any("load_md5_new_sequences$in_md5" in s for s in sql_calls)


def test_run_defaults_to_auto_release_type(monkeypatch):
    calls = _run_capture(monkeypatch)
    sql_calls = [sql for sql, _ in calls]

    # Release type is chosen per database from history: 'A' (auto), not a forced 'F'.
    assert ("SELECT rnc_update.prepare_releases(%s)", ("A",)) in calls
    assert "SELECT rnc_update.prepare_releases('F')" not in sql_calls

    # The per-database load still runs for each release returned by the pending query.
    assert ("SELECT rnc_update.new_update_release(%s, %s)", (9, 123)) in calls

    # do_checks runs exactly once, after the per-database loop -- but not last,
    # since _record_pending_stats runs after it.
    assert sum("do_checks(NULL::bigint)" in sql for sql in sql_calls) == 1
    assert sql_calls.index("SELECT rnc_load_xref.do_checks(NULL::bigint)") > max(
        i for i, s in enumerate(sql_calls) if "new_update_release" in s
    )


def test_run_force_full_forces_full_release_type(monkeypatch):
    calls = _run_capture(monkeypatch, force_full=True)

    # --force-full pins every database to a FULL release.
    assert ("SELECT rnc_update.prepare_releases(%s)", ("F",)) in calls
    assert ("SELECT rnc_update.prepare_releases(%s)", ("A",)) not in calls
