#!/usr/bin/env python3
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
Record / refresh the network cassettes used by the offline test suite.

Most of the time you do not need this: running any test with
``RNAC_TEST_ALLOW_NETWORK=1`` records any interaction that is not already
cached, so brand new network calls are captured automatically. This script is
for the cases that are *not* automatic:

* ``--refresh`` - delete the existing cassettes first so stale recordings are
  replaced with current data (record-on-miss never overwrites an existing
  cassette), and
* re-seeding the handful of fixtures that have to be maintained by hand because
  their upstream endpoint is gone (see ``MANUAL_SEEDS`` below).

Usage:
    uv run python tests/update_cassettes.py            # fill any missing cassettes
    uv run python tests/update_cassettes.py --refresh  # wipe and re-record everything
    uv run python tests/update_cassettes.py tests/databases/ols  # only these targets
"""

import argparse
import base64
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CASSETTE_DIR = ROOT / "tests" / "data" / "cassettes"
PHYLOGENY_DIR = ROOT / "tests" / "data" / "phylogeny"

# Test targets that talk to external services. Anything not listed still records
# correctly if you point the script at it explicitly, but these are the paths
# whose cassettes we maintain.
DEFAULT_TARGETS = [
    "tests/databases/pdb",
    "tests/databases/ncbi/gene",
    "tests/databases/gtrnadb/urls_test.py",
    "tests/databases/pirbase/fetch_test.py",
    "tests/databases/ena/dr_test.py",
    "tests/databases/europepmc/xml_test.py",
    "tests/databases/data/references_test.py",
    "tests/databases/ensembl/genomes/urls_test.py",
    "tests/databases/ensembl/metadata/databases_test.py",
    "tests/databases/ensembl/metadata/karyotypes_test.py",
    "tests/databases/sequence_ontology/tree_test.py",
    "tests/databases/data/rna_type_test.py",
    "tests/cli/ols_test.py",
    "tests/cli/europepmc_test.py",
    "tests/cli/pdb_test.py",
]

# Cassettes that cannot be recorded live because their upstream is gone, but for
# which we have a real historical payload vendored in the repo. These are
# (re)written after recording so that a --refresh does not leave them broken.
# Maps the cassette key (filename without extension) -> a seeder callable.
MANUAL_SEEDS = {
    # http://trna.ucsc.edu/.../export2019/ now 301s to https which 403s; the real
    # 2019 directory listing lives in data/gtrnadb/files.html.
    "requests/913de0502f52e7438e0d0fae718645e6b79cce7106a97546a65b0b859b159384": (
        lambda: _html_record(
            "https://trna.ucsc.edu/download/RNAcentral/export2019/",
            ROOT / "data" / "gtrnadb" / "files.html",
        )
    ),
}


def _html_record(url, html_path):
    return {
        "status_code": 200,
        "url": url,
        "reason": "OK",
        "encoding": "UTF-8",
        "headers": {"Content-Type": "text/html;charset=UTF-8"},
        "content_b64": base64.b64encode(html_path.read_bytes()).decode("ascii"),
    }


def _apply_manual_seeds():
    import json

    for key, builder in MANUAL_SEEDS.items():
        path = CASSETTE_DIR / f"{key}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(builder(), indent=2, sort_keys=True))
        print(f"  seeded manual cassette {key}")


def main():
    parser = argparse.ArgumentParser(
        description="Record or refresh the network cassettes used by the offline "
        "test suite. With no arguments it fills any missing cassettes; --refresh "
        "wipes and re-records them. See the module docstring for details."
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="delete existing cassettes (and the phylogeny cache) before recording",
    )
    parser.add_argument(
        "targets",
        nargs="*",
        help="optional pytest targets to record (defaults to the maintained set)",
    )
    args = parser.parse_args()

    targets = args.targets or DEFAULT_TARGETS

    if args.refresh:
        print("Removing existing cassettes and phylogeny cache...")
        shutil.rmtree(CASSETTE_DIR, ignore_errors=True)
        shutil.rmtree(PHYLOGENY_DIR, ignore_errors=True)

    env = dict(os.environ, RNAC_TEST_ALLOW_NETWORK="1")
    cmd = [
        sys.executable,
        "-m",
        "pytest",
        "-p",
        "no:cacheprovider",
        "-o",
        "addopts=",
        "-q",
        "--tb=no",
        *targets,
    ]
    print("Recording cassettes (live network)...")
    print("  " + " ".join(cmd))
    # Test assertions may still fail where upstream data has drifted; that does
    # not stop the cassettes from being recorded, so we ignore the exit code.
    subprocess.run(cmd, cwd=ROOT, env=env)

    print("Re-applying manual seeds...")
    _apply_manual_seeds()

    print("\nDone. Review changes with `git status tests/data` and re-run the")
    print("suite offline (no env var) to confirm there are no cassette misses.")


if __name__ == "__main__":
    main()
