# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

# Shared decoder for the active-sequence dump (files/ftp-export/sequences/active.sql).
# The dump is `COPY (...) TO STDOUT`, whose text format doubles backslashes, so a
# strict json.loads sometimes fails. Both the FASTA export (bin/json2fasta.py) and
# the parquet export (parquet.py) parse the dump through here so they can never drift.

import json
import logging

LOGGER = logging.getLogger(__name__)


def decode_entry(raw):
    try:
        return json.loads(raw)
    except json.JSONDecodeError as error:
        try:
            return json.loads(raw.replace("\\\\", "\\"))
        except json.JSONDecodeError:
            raise error


def parse(handle):
    """
    Yield one dict per line. The dump is `COPY (json_build_object(...)) TO STDOUT`,
    which emits exactly one compact JSON object per physical line (embedded newlines
    are escaped by COPY), so each line decodes on its own. We must NOT buffer across
    lines: a single undecodable line would otherwise accumulate the entire rest of
    the (multi-hundred-GB) dump into one string and OOM. Fail fast on a bad line
    instead, naming it, so the dump can be fixed rather than silently truncated.
    """
    for n, line in enumerate(handle, start=1):
        if not line.strip():
            continue
        try:
            yield decode_entry(line)
        except json.JSONDecodeError:
            LOGGER.error("Undecodable JSON on line %d: %.200r", n, line)
            raise


if __name__ == "__main__":
    import io

    # ponytail: covers the strict line and the COPY backslash-doubling fallback.
    entries = list(parse(io.StringIO('{"id": "a", "sequence": "ACGT"}\n')))
    assert entries == [{"id": "a", "sequence": "ACGT"}], entries
    doubled = list(parse(io.StringIO('{"id": "b", "description": "x\\\\ny", "sequence": "AC"}\n')))
    assert doubled[0]["description"] == "x\\ny", doubled

    # per-line streaming: each valid line yields immediately, no cross-line buffering
    good = "".join('{"id": "u%d"}\n' % i for i in range(3))
    assert [e["id"] for e in parse(io.StringIO(good))] == ["u0", "u1", "u2"]

    # a bad line must fail fast (naming the line), not silently swallow the rest
    try:
        list(parse(io.StringIO('{"id": "u0"}\n{bad\n{"id": "u2"}\n')))
        raise AssertionError("expected JSONDecodeError on the bad line")
    except json.JSONDecodeError:
        pass
    print("ok")
