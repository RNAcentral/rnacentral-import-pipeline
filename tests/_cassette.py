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
            key = _key("mysql", host, db, query, args)
            record = _load("mysql", key)
            if record is None:
                if not _ALLOW_NETWORK:
                    _miss("MySQL", f"{host}/{db}: {query}")
                cursor = self._connection._real.cursor()
                cursor.execute(query, args)
                rows = cursor.fetchall()
                cursor.close()
                record = {"rows": [list(r) for r in rows]}
                _save("mysql", key, record)
            self._rows = [tuple(r) for r in record["rows"]]
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

    class Row(list):
        """List-like, but also supports column-name access (like DictRow) and
        dict(row) (dict() treats anything with .keys() as a mapping)."""

        def __init__(self, values, columns):
            super().__init__(values)
            self._columns = columns or []

        def keys(self):
            return list(self._columns)

        def __getitem__(self, key):
            if isinstance(key, str):
                return list.__getitem__(self, self._columns.index(key))
            return list.__getitem__(self, key)

    class CassetteCursor:
        def __init__(self, connection):
            self._connection = connection
            self._rows = []
            self._columns = None
            self._pos = 0
            self.rowcount = -1

        def execute(self, query, params=None):
            # No DSN in the key: replay must work with whatever (or no)
            # PGDATABASE the current environment has, not just the one used
            # to record - there is only one real target in practice.
            key = _key("psycopg2", query, params)
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
            self._rows = [
                Row([_decode_cell(v) for v in row], self._columns)
                for row in record["rows"]
            ]
            self._pos = 0
            self.rowcount = record["rowcount"]
            return self.rowcount

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
