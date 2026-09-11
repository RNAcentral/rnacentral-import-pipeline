# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

"""
Record/replay cassettes for the external services the test suite talks to.

Each transport (HTTP via ``requests``, NCBI Entrez via urllib, FTP via
``ftplib``, Ensembl MySQL via ``pymysql``, and our own Postgres via
``psycopg2``/the ``psql`` CLI) is patched so that:

* in replay mode (the default) responses are served from vendored fixtures
  under ``tests/data/cassettes/`` and any un-recorded interaction raises, and
* in record mode (``RNAC_TEST_ALLOW_NETWORK=1``) the real call is made and the
  response written back to the cassette.

Because the patches sit at the transport boundary the real parsing code runs
against real (recorded) payloads - the fixtures are genuine API responses, just
captured once rather than fetched on every run.
"""

import base64
import datetime as dt
import decimal
import hashlib
import io
import json
import re
from pathlib import Path
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

# Hosts a Postgres cassette never wraps: tests/sql/conftest.py spins up its own
# throwaway local Postgres per run, so recording/replaying it would just pin
# stale schema-creation SQL against a database that no longer exists.
_LOCAL_DB_HOSTS = {"", "127.0.0.1", "localhost", "::1"}

CASSETTE_DIR = Path(__file__).parent.parent / "test-data" / "cassettes"

# Populated by install(); restore callables are run by uninstall().
_RESTORE = []
_ALLOW_NETWORK = False


class CassetteMiss(RuntimeError):
    """Raised when an interaction is not recorded and network is disabled."""


def _key(*parts):
    digest = hashlib.sha256()
    for part in parts:
        if part is None:
            part = b""
        elif isinstance(part, str):
            part = part.encode("utf-8")
        elif not isinstance(part, (bytes, bytearray)):
            part = str(part).encode("utf-8")
        digest.update(part)
        digest.update(b"\x1f")
    return digest.hexdigest()


def _path(kind, key):
    return CASSETTE_DIR / kind / f"{key}.json"


def _load(kind, key):
    path = _path(kind, key)
    if path.exists():
        return json.loads(path.read_text())
    return None


def _save(kind, key, record):
    path = _path(kind, key)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True))


def _miss(kind, description):
    raise CassetteMiss(
        f"No recorded {kind} interaction for {description}. Record it by running "
        "with RNAC_TEST_ALLOW_NETWORK=1."
    )


def _canonical_tokens(value):
    """Sort comma-separated tokens so set-ordered fields key consistently."""
    if "," in value:
        return ",".join(sorted(value.split(",")))
    return value


def _canonical_url(url):
    """
    Canonicalise a URL for use as a cassette key.

    Several callers build query parameters (and comma-joined value lists such as
    PDBe's ``fl=`` field list) from Python sets, whose iteration order varies
    between processes. Sorting the parameters and their comma-separated tokens
    keeps the record and replay keys stable.
    """
    parts = urlsplit(url)
    if not parts.query:
        return url
    params = [(k, _canonical_tokens(v)) for k, v in parse_qsl(parts.query)]
    params.sort()
    return urlunsplit(
        (parts.scheme, parts.netloc, parts.path, urlencode(params), parts.fragment)
    )


def _canonical_body(body):
    if body is None:
        return b""
    if isinstance(body, str):
        body = body.encode("utf-8")
    try:
        text = body.decode("utf-8")
    except (UnicodeDecodeError, AttributeError):
        return body
    # Comma-joined identifier lists (e.g. PDBe publication POSTs) are also
    # built from sets, so canonicalise the same way.
    if text and all(c.isalnum() or c in ",._-:" for c in text):
        return _canonical_tokens(text).encode("utf-8")
    return body


def _canonical_params(value):
    """
    Sort a homogeneous list/tuple of scalars for use in a cassette key.

    Several SQL callers build an id list for `= ANY(%s)`/`IN %s` from a
    Python set (e.g. a set comprehension over parsed rows), whose iteration
    order is hash-randomized per process - recording in one process and
    replaying in another would otherwise almost never hit the same key. Only
    a *homogeneous* collection is sorted (every element the same type): a
    positional args tuple like (min_id, max_id) must keep its order, and
    mixed-type tuples are never how this codebase builds an id list.
    """
    if isinstance(value, (list, tuple)):
        items = [_canonical_params(v) for v in value]
        if (
            len(items) > 1
            and len({type(v) for v in items}) == 1
            and isinstance(items[0], (str, int, float, bool))
        ):
            return sorted(items)
        return items
    if isinstance(value, dict):
        return {k: _canonical_params(v) for k, v in value.items()}
    return value


class Row(list):
    """
    List-like, but also supports column-name access (like DictRow/DictCursor)
    and dict(row) (dict() treats anything with .keys() as a mapping).

    Module-level rather than nested in an _install_* function: production
    code sometimes pickles a row (e.g. rnacentral_pipeline.rnacentral.lookup's
    write_mapping), and pickle needs a class importable by its qualified name
    - a class local to a function isn't.
    """

    def __init__(self, values, columns):
        super().__init__(values)
        self._columns = columns or []

    def keys(self):
        return list(self._columns)

    def __getitem__(self, key):
        if isinstance(key, str):
            return list.__getitem__(self, self._columns.index(key))
        return list.__getitem__(self, key)

    def __setitem__(self, key, value):
        if isinstance(key, str):
            return list.__setitem__(self, self._columns.index(key), value)
        return list.__setitem__(self, key, value)


# ---------------------------------------------------------------------------
# requests (HTTP/HTTPS)
# ---------------------------------------------------------------------------
def _install_requests():
    import requests.exceptions as rexc
    from requests.adapters import HTTPAdapter
    from requests.models import Response
    from requests.structures import CaseInsensitiveDict

    real_send = HTTPAdapter.send

    def send(self, request, **kwargs):
        key = _key(
            "requests",
            request.method,
            _canonical_url(request.url),
            _canonical_body(request.body),
        )
        record = _load("requests", key)
        if record is None:
            if not _ALLOW_NETWORK:
                _miss("HTTP", f"{request.method} {request.url}")
            try:
                resp = real_send(self, request, **kwargs)
            except rexc.RequestException as err:
                record = {"error": type(err).__name__, "message": str(err)}
                _save("requests", key, record)
            else:
                content = resp.content
                record = {
                    "status_code": resp.status_code,
                    "url": resp.url,
                    "reason": resp.reason or "",
                    "encoding": resp.encoding,
                    # Store all headers (not just Content-Type): requests needs
                    # Location to follow redirects on replay, and the redirect
                    # target is recorded as a separate interaction.
                    "headers": dict(resp.headers),
                    "content_b64": base64.b64encode(content).decode("ascii"),
                }
                _save("requests", key, record)

        if "error" in record:
            exc_cls = getattr(rexc, record["error"], rexc.RequestException)
            raise exc_cls(record["message"])

        response = Response()
        response.status_code = record["status_code"]
        response.url = record["url"]
        response.reason = record["reason"]
        response.encoding = record["encoding"]
        response.headers = CaseInsensitiveDict(record["headers"])
        response._content = base64.b64decode(record["content_b64"])
        response.request = request
        return response

    HTTPAdapter.send = send
    _RESTORE.append(lambda: setattr(HTTPAdapter, "send", real_send))


# ---------------------------------------------------------------------------
# Bio.Entrez (NCBI E-utilities over urllib)
# ---------------------------------------------------------------------------
def _install_entrez():
    try:
        from Bio import Entrez
    except Exception:
        return

    real_open = Entrez._open

    def fake_open(request, *args, **kwargs):
        url = getattr(request, "full_url", None) or request.get_full_url()
        data = getattr(request, "data", None)
        key = _key("entrez", url, data)
        record = _load("entrez", key)
        if record is None:
            if not _ALLOW_NETWORK:
                _miss("Entrez", url)
            handle = real_open(request, *args, **kwargs)
            # _open hands back a binary handle for XML payloads and wraps plain
            # text responses in a TextIOWrapper; mirror that on replay so that
            # both Entrez.read (binary) and SeqIO.parse (text) stay happy.
            is_text = isinstance(handle, io.TextIOWrapper)
            raw = handle.read()
            handle.close()
            content = raw.encode("utf-8") if isinstance(raw, str) else raw
            record = {
                "url": url,
                "is_text": is_text,
                "content_b64": base64.b64encode(content).decode("ascii"),
            }
            _save("entrez", key, record)
        stream = io.BytesIO(base64.b64decode(record["content_b64"]))
        if record.get("is_text"):
            return io.TextIOWrapper(stream, encoding="utf-8")
        return stream

    # Entrez._open is decorated with @function_with_previous, which records the
    # timestamp of the last call on the function object itself and references it
    # through the module global. Replacing that global means the original (which
    # we still call in record mode) reads ``.previous`` off our replacement, so
    # it must carry the attribute too.
    fake_open.previous = getattr(real_open, "previous", 0.0)
    Entrez._open = fake_open
    _RESTORE.append(lambda: setattr(Entrez, "_open", real_open))


# ---------------------------------------------------------------------------
# ftplib (Ensembl/NCBI FTP file retrieval)
# ---------------------------------------------------------------------------
def _install_ftp():
    import ftplib

    real_ftp = ftplib.FTP

    class CassetteFTP:
        def __init__(self, host="", *args, **kwargs):
            self._host = host
            self._cwd = "/"
            self._real = None
            if _ALLOW_NETWORK:
                self._real = real_ftp(host, *args, **kwargs)

        # Context manager -------------------------------------------------
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.close()
            return False

        # Navigation ------------------------------------------------------
        def login(self, *args, **kwargs):
            if self._real is not None:
                return self._real.login(*args, **kwargs)
            return "230 Login successful."

        def cwd(self, path):
            if self._real is not None:
                result = self._real.cwd(path)
            else:
                result = "250 OK."
            if path.startswith("/"):
                self._cwd = path
            else:
                self._cwd = (self._cwd.rstrip("/") + "/" + path).replace("//", "/")
            return result

        def quit(self):
            return self.close()

        def close(self):
            if self._real is not None:
                try:
                    self._real.quit()
                except Exception:
                    self._real.close()
                self._real = None

        # Retrieval -------------------------------------------------------
        def retrlines(self, cmd, callback=None):
            key = _key("ftp", self._host, self._cwd, cmd)
            record = _load("ftp", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("FTP", f"{self._host}:{self._cwd} {cmd}")
                lines = []
                try:
                    self._real.retrlines(cmd, lines.append)
                    record = {"lines": lines}
                except ftplib.Error as err:
                    record = {"error": type(err).__name__, "message": str(err)}
                _save("ftp", key, record)
            if "error" in record:
                raise getattr(ftplib, record["error"])(record["message"])
            for line in record["lines"]:
                if callback is not None:
                    callback(line)
            return "226 Transfer complete."

        def retrbinary(self, cmd, callback, blocksize=8192, *args, **kwargs):
            key = _key("ftp", self._host, self._cwd, cmd)
            record = _load("ftp", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("FTP", f"{self._host}:{self._cwd} {cmd}")
                chunks = []
                try:
                    self._real.retrbinary(
                        cmd, chunks.append, blocksize, *args, **kwargs
                    )
                    payload = b"".join(chunks)
                    record = {"content_b64": base64.b64encode(payload).decode("ascii")}
                except ftplib.Error as err:
                    record = {"error": type(err).__name__, "message": str(err)}
                _save("ftp", key, record)
            if "error" in record:
                raise getattr(ftplib, record["error"])(record["message"])
            callback(base64.b64decode(record["content_b64"]))
            return "226 Transfer complete."

        def nlst(self, *args):
            cmd = "NLST " + " ".join(args)
            key = _key("ftp", self._host, self._cwd, cmd)
            record = _load("ftp", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("FTP", f"{self._host}:{self._cwd} {cmd}")
                record = {"names": self._real.nlst(*args)}
                _save("ftp", key, record)
            return record["names"]

    ftplib.FTP = CassetteFTP
    _RESTORE.append(lambda: setattr(ftplib, "FTP", real_ftp))


# ---------------------------------------------------------------------------
# pymysql (Ensembl public MySQL)
# ---------------------------------------------------------------------------
def _install_pymysql():
    try:
        import pymysql
    except Exception:
        return

    real_connection = pymysql.connections.Connection

    class CassetteCursor:
        def __init__(self, connection):
            self._connection = connection
            self._rows = []

        def execute(self, query, args=None):
            host = self._connection._cassette_host
            db = self._connection._cassette_db
            key = _key("mysql", host, db, query, _canonical_params(args))
            record = _load("mysql", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("MySQL", f"{host}/{db}: {query}")
                cursor = self._connection._real.cursor()
                cursor.execute(query, args)
                columns = (
                    [d[0] for d in cursor.description] if cursor.description else None
                )
                rows = cursor.fetchall()
                cursor.close()
                record = {"columns": columns, "rows": [list(r) for r in rows]}
                _save("mysql", key, record)
            columns = record.get("columns")
            self._rows = [Row(r, columns) for r in record["rows"]]
            return len(self._rows)

        def fetchall(self):
            return list(self._rows)

        def fetchone(self):
            return self._rows.pop(0) if self._rows else None

        def close(self):
            pass

        def __iter__(self):
            return iter(self._rows)

    class CassetteConnection:
        def __init__(self, *args, **kwargs):
            self._cassette_host = kwargs.get("host", "")
            self._cassette_db = kwargs.get("database") or kwargs.get("db")
            self._real = None
            if _ALLOW_NETWORK:
                self._real = real_connection(*args, **kwargs)

        def cursor(self, cursor=None):
            return CassetteCursor(self)

        def select_db(self, db):
            self._cassette_db = db
            if self._real is not None:
                self._real.select_db(db)

        def close(self):
            if self._real is not None:
                self._real.close()
                self._real = None

    pymysql.connections.Connection = CassetteConnection
    pymysql.Connection = CassetteConnection
    pymysql.connect = CassetteConnection
    _RESTORE.append(
        lambda: (
            setattr(pymysql.connections, "Connection", real_connection),
            setattr(pymysql, "Connection", real_connection),
            setattr(pymysql, "connect", real_connection),
        )
    )


# ---------------------------------------------------------------------------
# psycopg2 (our own Postgres, e.g. rnacentral_pipeline.rnacentral.lookup)
# ---------------------------------------------------------------------------
def _encode_cell(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return {"__array__": [_encode_cell(v) for v in value]}
    if isinstance(value, dict):
        return {"__dict__": {k: _encode_cell(v) for k, v in value.items()}}
    if isinstance(value, decimal.Decimal):
        return {"__decimal__": str(value)}
    if isinstance(value, (dt.datetime, dt.date)):
        return {"__datetime__": value.isoformat()}
    if isinstance(value, (bytes, bytearray, memoryview)):
        return {"__bytes__": base64.b64encode(bytes(value)).decode("ascii")}
    return {"__repr__": str(value)}


def _decode_cell(value):
    if isinstance(value, dict):
        if "__array__" in value:
            return [_decode_cell(v) for v in value["__array__"]]
        if "__dict__" in value:
            return {k: _decode_cell(v) for k, v in value["__dict__"].items()}
        if "__decimal__" in value:
            return decimal.Decimal(value["__decimal__"])
        if "__datetime__" in value:
            raw = value["__datetime__"]
            return (
                dt.datetime.fromisoformat(raw)
                if "T" in raw
                else dt.date.fromisoformat(raw)
            )
        if "__bytes__" in value:
            return base64.b64decode(value["__bytes__"])
        if "__repr__" in value:
            return value["__repr__"]
    return value


def _install_psycopg2():
    try:
        import psycopg2
        import psycopg2.errors
        import psycopg2.extensions
    except Exception:
        return

    real_connect = psycopg2.connect

    def _host_of_dsn(args, kwargs):
        try:
            if args:
                parsed = psycopg2.extensions.parse_dsn(args[0])
            else:
                parsed = kwargs
            return parsed.get("host", "")
        except Exception:
            return ""

    class CassetteCursor:
        def __init__(self, connection):
            self._connection = connection
            self._rows = []
            self._columns = None
            self._pos = 0
            self.rowcount = -1
            # Real psycopg2 cursors always expose .description (a sequence of
            # 7-tuples, one per column - only the first element, the name, is
            # ever populated here) after execute(); DB-API consumers that
            # don't go through fetchall() (e.g. polars' row-wise fallback)
            # read this directly rather than asking the cursor for columns.
            self.description = None

        def execute(self, query, params=None, *, vars=None, parameters=None):
            # Real psycopg2's execute() is a C method, so inspect.signature()
            # raises on it - callers that branch on introspection (e.g.
            # polars.read_database, which decides positional-vs-keyword args
            # this way) take a *different* branch against this pure-Python
            # stand-in and end up calling with vars=/parameters= instead of
            # the positional/params call real psycopg2 would receive. Accept
            # all three spellings so either branch works the same.
            if params is None:
                params = vars if vars is not None else parameters
            # execute_batch/execute_values (psycopg2.extras) mogrify() each
            # row and pass the joined result to execute() as bytes, unlike
            # every other caller's plain str query.
            query_text = query.decode("utf-8") if isinstance(query, bytes) else query
            # Writes/DDL (CREATE TEMP TABLE, INSERT, SET, ...) are never
            # cached: their point is a session-local side effect, not a
            # reusable result, and a fixture that recreates the same temp
            # table per-test would otherwise get a cache hit on test 2+ that
            # skips the real CREATE - leaving that test's session without the
            # table it just "created", silently falling through to whatever
            # real table of the same name exists. Always run for real while
            # recording; a no-op offline, since replay has no live session to
            # affect and the point was only ever to enable later real reads.
            if not re.match(r"^\s*(SELECT|WITH)\b", query_text, re.IGNORECASE):
                self._columns = None
                self._rows = []
                self.description = None
                self.rowcount = 0
                if self._connection._real is not None:
                    real_cur = self._connection._real.cursor()
                    real_cur.execute(query, params)
                    self.rowcount = real_cur.rowcount
                    real_cur.close()
                # Real psycopg2 execute() always returns None (results are
                # read back off the cursor, not off a return value) - a
                # caller that branches on the return value (e.g. polars.
                # read_database's "did this execute in-place" check) needs
                # that exact contract, not our rowcount.
                return None

            # No DSN in the key: replay must work with whatever (or no)
            # PGDATABASE the current environment has, not just the one used
            # to record - there is only one real target in practice.
            key = _key("psycopg2", query, _canonical_params(params))
            record = _load("psycopg2", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("Postgres", query)
                real_cur = self._connection._real.cursor()
                # query/params are stored too, alongside the result, purely so
                # a drift-checker can walk the cassette dir and re-ask the same
                # question later - replay itself never reads these two keys.
                base = {"query": query, "params": _encode_cell(params)}
                try:
                    real_cur.execute(query, params)
                except psycopg2.Error as err:
                    record = {
                        **base,
                        "error": type(err).__name__,
                        "message": str(err),
                    }
                else:
                    columns = (
                        [d[0] for d in real_cur.description]
                        if real_cur.description
                        else None
                    )
                    rows = real_cur.fetchall() if real_cur.description else []
                    record = {
                        **base,
                        "columns": columns,
                        "rowcount": real_cur.rowcount,
                        "rows": [[_encode_cell(v) for v in row] for row in rows],
                    }
                real_cur.close()
                _save("psycopg2", key, record)

            if "error" in record:
                exc_cls = getattr(psycopg2.errors, record["error"], None) or getattr(
                    psycopg2, record["error"], psycopg2.Error
                )
                raise exc_cls(record["message"])

            self._columns = record["columns"]
            self.description = (
                [(c, None, None, None, None, None, None) for c in self._columns]
                if self._columns
                else None
            )
            self._rows = [
                Row([_decode_cell(v) for v in row], self._columns)
                for row in record["rows"]
            ]
            self._pos = 0
            self.rowcount = record["rowcount"]
            return None

        def fetchall(self):
            rows = self._rows[self._pos :]
            self._pos = len(self._rows)
            return rows

        def fetchone(self):
            if self._pos >= len(self._rows):
                return None
            row = self._rows[self._pos]
            self._pos += 1
            return row

        def copy_expert(self, sql, file, size=8192):
            # COPY ... TO STDOUT streams raw text into a file-like object -
            # a completely different shape of call than execute()/fetchall(),
            # so it gets its own cache entry rather than reusing that path.
            key = _key("psycopg2", "COPY", sql)
            record = _load("psycopg2", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("Postgres", sql)
                real_cur = self._connection._real.cursor()
                buf = io.StringIO()
                real_cur.copy_expert(sql, buf, size)
                real_cur.close()
                record = {"query": sql, "data": buf.getvalue()}
                _save("psycopg2", key, record)
            file.write(record["data"])

        def mogrify(self, query, vars=None):
            # Used by psycopg2.extras.execute_batch/execute_values to build a
            # combined multi-row statement before a single execute() call.
            # Never cached: like other writes, the point is safely interpolated
            # SQL text to actually run while recording, not a reusable result -
            # offline there is nothing to run it against, so the exact bytes
            # don't matter (the caller's follow-up execute() is itself a no-op
            # without a live connection).
            if self._connection._real is not None:
                real_cur = self._connection._real.cursor()
                try:
                    return real_cur.mogrify(query, vars)
                finally:
                    real_cur.close()
            return b""

        def close(self):
            pass

        def __iter__(self):
            return iter(self._rows)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.close()
            return False

    class CassetteConnection:
        def __init__(self, *args, **kwargs):
            self._real = None
            if _ALLOW_NETWORK:
                self._real = real_connect(*args, **kwargs)

        def cursor(self, *args, **kwargs):
            return CassetteCursor(self)

        def commit(self):
            if self._real is not None:
                self._real.commit()

        def rollback(self):
            if self._real is not None:
                self._real.rollback()

        def set_session(self, *args, **kwargs):
            if self._real is not None:
                self._real.set_session(*args, **kwargs)

        def set_isolation_level(self, level):
            if self._real is not None:
                self._real.set_isolation_level(level)

        @property
        def autocommit(self):
            return self._real.autocommit if self._real is not None else False

        @autocommit.setter
        def autocommit(self, value):
            if self._real is not None:
                self._real.autocommit = value

        def close(self):
            if self._real is not None:
                self._real.close()
                self._real = None

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.close()
            return False

    def connect(*args, **kwargs):
        if _host_of_dsn(args, kwargs) in _LOCAL_DB_HOSTS:
            return real_connect(*args, **kwargs)
        return CassetteConnection(*args, **kwargs)

    psycopg2.connect = connect
    _RESTORE.append(lambda: setattr(psycopg2, "connect", real_connect))


# ---------------------------------------------------------------------------
# psql CLI (tests/helpers.py shells out to it rather than using psycopg2)
# ---------------------------------------------------------------------------
def _psql_target_is_local(db_arg):
    """The last positional arg to `psql` is a dbname or a full connection URI."""
    try:
        import psycopg2.extensions

        return psycopg2.extensions.parse_dsn(db_arg).get("host", "") in _LOCAL_DB_HOSTS
    except Exception:
        # A bare dbname (no `://`) connects via the local default host/socket.
        return "://" not in db_arg


def _install_psql_subprocess():
    import subprocess

    real_run = subprocess.run

    def fake_run(args, *a, **kwargs):
        argv = list(args) if isinstance(args, (list, tuple)) else None
        if (
            argv
            and len(argv) >= 2
            and Path(str(argv[0])).name == "psql"
            and "-f" in argv
            and not _psql_target_is_local(argv[-1])
        ):
            sql_path = Path(argv[argv.index("-f") + 1])
            sql_text = sql_path.read_text() if sql_path.exists() else ""
            key = _key("psql", sql_text)
            record = _load("psql", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("Postgres", f"psql -f {sql_path.name}")
                result = real_run(args, *a, **kwargs)
                stdout = result.stdout
                if isinstance(stdout, bytes):
                    stdout = stdout.decode("utf-8", "replace")
                # sql is stored alongside the result purely so a drift-checker
                # can re-run the same query later - replay never reads it.
                record = {
                    "sql": sql_text,
                    "returncode": result.returncode,
                    "stdout": stdout or "",
                }
                _save("psql", key, record)
            return subprocess.CompletedProcess(
                args, record["returncode"], stdout=record["stdout"]
            )
        return real_run(args, *a, **kwargs)

    subprocess.run = fake_run
    _RESTORE.append(lambda: setattr(subprocess, "run", real_run))


def install(allow_network):
    global _ALLOW_NETWORK
    _ALLOW_NETWORK = allow_network
    _install_requests()
    _install_entrez()
    _install_ftp()
    _install_pymysql()
    _install_psycopg2()
    _install_psql_subprocess()


def uninstall():
    while _RESTORE:
        restore = _RESTORE.pop()
        try:
            restore()
        except Exception:
            pass
