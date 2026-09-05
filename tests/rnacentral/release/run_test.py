import psycopg2.extensions

from rnacentral_pipeline.rnacentral.release import run


class FakeCursor:
    def __init__(self, conn):
        self.conn = conn
        self.calls = []
        self.opts_calls = []
        self._last = ""

    def execute(self, sql, params=None):
        self.calls.append((sql, params))
        # Each step opens its own connection, so the options recorded at connect
        # time are the budget this statement runs under.
        self.opts_calls.append((sql, self.conn.current_options))
        self._last = sql
        # Stand in for a RAISE NOTICE arriving from the running function.
        if sql == "SELECT rnc_update.update_rnc_accessions()":
            self.conn.notices.append("NOTICE:  Upserting rn [1, 15000000)\n")

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
        self.cursor_obj = FakeCursor(self)
        self.commits = 0
        self.autocommit = False
        self.notices = []
        self.connect_kwargs = []
        self.current_options = None
        self.closed = 0

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

    def poll(self):
        return psycopg2.extensions.POLL_OK

    def fileno(self):
        raise AssertionError("select() reached even though poll() said POLL_OK")

    def close(self):
        self.closed += 1


def _run_release(monkeypatch, **kwargs):
    conn = FakeConnection()

    # _connect passes keepalive/options kwargs, so the stub must accept them.
    def fake_connect(*args, **opts):
        conn.connect_kwargs.append(opts)
        conn.current_options = opts.get("options")
        return conn

    monkeypatch.setattr(run.psycopg2, "connect", fake_connect)
    run.run("postgres://example", **kwargs)
    return conn


def _run_capture(monkeypatch, **kwargs):
    return _run_release(monkeypatch, **kwargs).cursor_obj.calls


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


def test_run_logs_progress_for_every_step(monkeypatch, caplog):
    """A release run that logs nothing is indistinguishable from a wedged one."""
    with caplog.at_level("INFO", logger=run.LOGGER.name):
        _run_capture(monkeypatch)

    messages = [r.getMessage() for r in caplog.records]

    # Each long-running step brackets itself, so the last START with no matching
    # DONE is the step currently in progress.
    for step in (
        "apply_functions",
        "update_rnc_accessions",
        "update_literature_references",
        "prepare_releases",
        "new_update_release(dbid=9, rid=123) [1/1]",
        "do_checks (once, post-loop)",
    ):
        assert f"START  {step}" in messages, step
        assert any(m.startswith(f"DONE   {step} (") for m in messages), step

    assert "DONE   release run" in messages[-1]


def test_run_streams_server_notices_into_the_log(monkeypatch, caplog):
    """
    The RAISE NOTICE lines inside the release functions are the only view into a
    step that runs for hours, and a synchronous execute() hides them until the
    step finishes.
    """
    with caplog.at_level("INFO", logger=run.LOGGER.name):
        conn = _run_release(monkeypatch)

    messages = [r.getMessage() for r in caplog.records]

    # Tagged with the step it came from, so an interleaved log stays readable.
    assert "update_rnc_accessions | NOTICE:  Upserting rn [1, 15000000)" in messages

    # Notices only materialise on a polled (async) connection, and every one of
    # them is drained rather than left to pile up.
    assert any(opts.get("async_") == 1 for opts in conn.connect_kwargs)
    assert conn.notices == []


def _options_for(conn, predicate):
    # Match what run() issues, not function-deployment SQL: those sources mention
    # the same table and index names, so a substring search finds them first.
    return next(opts for sql, opts in conn.cursor_obj.opts_calls if predicate(sql))


def test_fk4_validation_runs_on_the_high_mem_connection(monkeypatch):
    """It joins a whole xref partition against rna, so not on the default budget."""
    conn = _run_release(monkeypatch)

    default_opts = _options_for(
        conn, lambda sql: sql == "SELECT rnc_update.update_rnc_accessions()"
    )
    validate_opts = _options_for(conn, lambda sql: sql.startswith("ALTER TABLE xref_p"))

    assert validate_opts != default_opts
    assert "work_mem=256MB" in validate_opts
