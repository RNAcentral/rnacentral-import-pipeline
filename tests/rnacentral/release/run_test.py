from rnacentral_pipeline.rnacentral.release import run


class FakeCursor:
    def __init__(self):
        self.calls = []
        self._results = [(9, 123)]

    def execute(self, sql, params=None):
        self.calls.append((sql, params))

    def fetchall(self):
        return self._results

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


def test_run_patches_functions_and_checks_once(monkeypatch):
    conn = FakeConnection()
    # _connect passes keepalive/options kwargs, so the stub must accept them.
    monkeypatch.setattr(run.psycopg2, "connect", lambda *a, **k: conn)

    run.run("postgres://example")

    sql_calls = [sql for sql, _ in conn.cursor_obj.calls]

    # do_pel_exchange is patched first for declarative partitioning support.
    assert "CREATE OR REPLACE FUNCTION rnc_load_xref.do_pel_exchange" in sql_calls[0]

    # load_xref is patched to drop the per-database do_checks call.
    assert "CREATE OR REPLACE FUNCTION rnc_load_xref.load_xref" in sql_calls[1]
    assert "perform rnc_load_xref.do_checks" not in sql_calls[1]

    # The per-database load still runs for each release returned by TO_RELEASE.
    assert ("SELECT rnc_update.new_update_release(%s, %s)", (9, 123)) in conn.cursor_obj.calls

    # do_checks runs exactly once, after the per-database loop.
    assert sql_calls[-1] == "SELECT rnc_load_xref.do_checks(NULL::bigint)"
    assert sum("do_checks(NULL::bigint)" in sql for sql in sql_calls) == 1
