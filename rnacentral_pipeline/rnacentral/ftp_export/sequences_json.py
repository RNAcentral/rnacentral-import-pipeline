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
    """Yield one dict per JSON line, buffering across lines for the rare multi-line entry."""
    pending = ""
    for line in handle:
        pending += line
        try:
            yield decode_entry(pending)
            pending = ""
        except json.JSONDecodeError:
            LOGGER.warning("Could not decode JSON entry yet, buffering additional input")
            continue
    if pending.strip():
        yield decode_entry(pending)


if __name__ == "__main__":
    import io

    # ponytail: covers the strict line and the COPY backslash-doubling fallback.
    entries = list(parse(io.StringIO('{"id": "a", "sequence": "ACGT"}\n')))
    assert entries == [{"id": "a", "sequence": "ACGT"}], entries
    doubled = list(parse(io.StringIO('{"id": "b", "description": "x\\\\ny", "sequence": "AC"}\n')))
    assert doubled[0]["description"] == "x\\ny", doubled
    print("ok")
