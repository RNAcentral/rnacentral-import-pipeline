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

"""
Static checks on the QC report wiring, which otherwise only runs hours into a
release: a missing psql var fails there, and a missing database just goes quiet.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
QC_WORKFLOW = ROOT / "workflows" / "utils" / "qc.nf"
MAIN_WORKFLOW = ROOT / "main.nf"
QC_SQL_DIR = ROOT / "files" / "qc"

PROCESS_RE = re.compile(r"^process\s+(\w+)\s*\{", re.MULTILINE)
WORKFLOW_RE = re.compile(r"^workflow\s+(\w+)\s*\{", re.MULTILINE)
PSQL_VAR_RE = re.compile(r"-v\s+(\w+)=")
QC_SQL_RE = re.compile(r"file\('files/qc/([\w-]+\.sql)'\)")
# :'quoted' and :bare references, but not ::casts.
SQL_VAR_RE = re.compile(r"(?<!:):'?([a-z_][a-z0-9_]*)'?")
SQL_IF_RE = re.compile(r"^\s*\\if\s+:\{\?(\w+)\}")
SQL_SET_RE = re.compile(r"^\s*\\set\s+(\w+)")
ENSEMBL_DIVISION_RE = re.compile(r"(\w+):\s*'ensembl\w*'")
MAIN_DIVISION_LIST_RE = re.compile(r"\[((?:'\w+',?\s*)+)\]\.any\s*\{\s*d\s*->")


def blocks(text, pattern):
    """Split a Nextflow file into {name: body} on the given top-level pattern."""
    starts = [(m.group(1), m.start()) for m in pattern.finditer(text)]
    boundaries = sorted(
        m.start() for m in list(PROCESS_RE.finditer(text)) + list(WORKFLOW_RE.finditer(text))
    )
    found = {}
    for name, start in starts:
        after = [b for b in boundaries if b > start]
        found[name] = text[start : after[0] if after else len(text)]
    return found


@pytest.fixture(scope="module")
def qc_nf():
    return QC_WORKFLOW.read_text()


@pytest.fixture(scope="module")
def invocations(qc_nf):
    """
    Every (process, sql file) pair the QC workflows run, read off the wiring
    rather than hard-coded, so a new stage is covered without touching this.
    """
    processes = blocks(qc_nf, PROCESS_RE)
    pairs = []
    for _name, body in blocks(qc_nf, WORKFLOW_RE).items():
        sql_files = QC_SQL_RE.findall(body)
        called = [p for p in processes if re.search(rf"\b{p}\b", body)]
        for sql in sql_files:
            for process in called:
                pairs.append((process, sql))
    assert pairs
    return sorted(set(pairs))


def passed_vars(qc_nf, process):
    return set(PSQL_VAR_RE.findall(blocks(qc_nf, PROCESS_RE)[process]))


def required_vars(sql_text, defined):
    """
    Variables the SQL needs from the caller: skips `\\if :{?var}` branches that
    cannot run, and counts a var the file `\\set`s itself as already satisfied.
    """
    defined = set(defined)
    required = set()
    active = True
    depth = 0

    for line in sql_text.splitlines():
        guard = SQL_IF_RE.match(line)
        if guard:
            depth += 1
            assert depth == 1, "nested psql conditionals are not handled"
            active = guard.group(1) in defined
            continue
        if re.match(r"^\s*\\else\b", line):
            active = not active
            continue
        if re.match(r"^\s*\\endif\b", line):
            depth -= 1
            active = True
            continue
        if not active:
            continue

        assignment = SQL_SET_RE.match(line)
        if assignment:
            defined.add(assignment.group(1))
            continue

        required.update(SQL_VAR_RE.findall(line.split("--")[0]))

    return required - defined


def test_every_qc_sql_file_is_wired_into_a_workflow(invocations):
    on_disk = {path.name for path in QC_SQL_DIR.glob("*.sql")}
    wired = {sql for _process, sql in invocations}
    assert on_disk == wired


def test_wired_qc_sql_files_exist(invocations):
    missing = [sql for _process, sql in invocations if not (QC_SQL_DIR / sql).exists()]
    assert missing == []


def test_each_process_passes_every_variable_its_sql_requires(qc_nf, invocations):
    problems = []
    for process, sql in invocations:
        given = passed_vars(qc_nf, process)
        needed = required_vars((QC_SQL_DIR / sql).read_text(), given)
        if needed - given:
            problems.append(f"{process} runs {sql} without {sorted(needed - given)}")

    assert problems == []


def test_optional_variables_are_guarded_rather_than_required(qc_nf, invocations):
    """
    A var with a documented default must not become mandatory: the SQL has to
    keep working when the process does not pass it.
    """
    problems = []
    for process, sql in invocations:
        text = (QC_SQL_DIR / sql).read_text()
        guarded = set(SQL_IF_RE.findall("\n".join(text.splitlines())))
        for var in guarded:
            if required_vars(text, passed_vars(qc_nf, process) - {var}) & {var}:
                problems.append(f"{sql}: {var} is guarded but still required")

    assert problems == []


def test_analyze_snapshot_does_not_need_the_run_start_variable():
    """
    The before-snapshot passes only -v snapshot=1, so nothing outside that
    branch may demand a variable it does not have.
    """
    text = (QC_SQL_DIR / "analyze.sql").read_text()
    assert required_vars(text, {"snapshot"}) == set()
    # The reporting run takes the other branch, which does need run_start.
    assert required_vars(text, set()) == {"run_start"}


def test_qc_knows_every_ensembl_division_main_runs(qc_nf):
    """
    qc.nf keeps its own division -> rnc_database.descr map, so a division added
    to main.nf but not here goes missing from the import report without error.
    """
    qc_divisions = set(ENSEMBL_DIVISION_RE.findall(qc_nf))
    main_lists = MAIN_DIVISION_LIST_RE.findall(MAIN_WORKFLOW.read_text())
    assert main_lists, "could not find the ensembl division list in main.nf"

    for raw in main_lists:
        assert set(re.findall(r"'(\w+)'", raw)) == qc_divisions
