from rnacentral_pipeline.rnacentral.release import run


class FakeCursor:
    def __init__(self):
        self.calls = []
        self._results = [(9, 123)]

    def execute(self, sql, params=None):
        self.calls.append((sql, params))

    def fetchall(self):
        return self._results


class FakeConnection:
    def __init__(self):
        self.cursor_obj = FakeCursor()
        self.commits = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def cursor(self):
        return self.cursor_obj

    def commit(self):
        self.commits += 1


def test_run_patches_xref_exchange_before_release(monkeypatch):
    conn = FakeConnection()
    monkeypatch.setattr(run.psycopg2, "connect", lambda _: conn)

    run.run("postgres://example")

    sql_calls = [sql for sql, _ in conn.cursor_obj.calls]
    assert "CREATE OR REPLACE FUNCTION rnc_load_xref.do_pel_exchange" in sql_calls[1]
    assert sql_calls[-1] == "SELECT rnc_update.new_update_release(%s, %s)"
    assert conn.cursor_obj.calls[-1][1] == (9, 123)
    assert conn.commits == 1
