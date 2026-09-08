# -*- coding: utf-8 -*-

"""
``precompute:process_range`` was killed by SLURM 728 tasks into a 1844-task run
and took the whole pipeline with it. The kill reached nextflow as

    Process `precompute:process_range (1175001-1200001)` terminated for an
    unknown reason -- Likely it has been terminated by the external system
    Command exit status: -

i.e. with no exit status at all, which is what nextflow reports when the batch
system removes the job step before ``.exitcode`` is written. The
``task.exitStatus in 137..140`` idiom used elsewhere in this repo cannot fire on
that, so the retry has to key off ``task.attempt``, and the memory request has
to grow with it or a re-run just dies at the same size.

``range.memory`` is multiplied by ``task.attempt``, so every config that sets it
must use a nextflow memory literal - a plain string would be *repeated* by
groovy's ``String * int`` and yield garbage like ``8GB8GB8GB``.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRECOMPUTE = (ROOT / "precompute.nf").read_text()
CONFIGS = sorted(ROOT.glob("**/*.config"))

RANGE_MEMORY = re.compile(r"^\s*range\.memory\s*=\s*(.+?)\s*$", re.MULTILINE)


def process_range_directives() -> str:
    body = PRECOMPUTE.split("process process_range {", 1)[1]
    return body.split("input:", 1)[0]


def test_process_range_retries_without_an_exit_status():
    directives = process_range_directives()

    strategy = [
        line.strip()
        for line in directives.splitlines()
        if line.strip().startswith("errorStrategy")
    ]

    assert strategy, "process_range has no errorStrategy: one SLURM kill ends the run"
    assert "task.attempt" in strategy[0], (
        f"process_range errorStrategy does not key off task.attempt: {strategy[0]!r} - "
        "an external kill arrives with no exit status, so an exitStatus guard "
        "never retries it"
    )
    assert "exitStatus" not in strategy[0], (
        f"process_range errorStrategy guards on exitStatus: {strategy[0]!r} - "
        "the failure it must survive reports no exit status at all"
    )


def test_process_range_memory_grows_with_the_attempt():
    directives = process_range_directives()

    memory = [
        line.strip()
        for line in directives.splitlines()
        if line.strip().startswith("memory")
    ]

    assert memory, "process_range has no memory directive"
    assert "task.attempt" in memory[0], (
        f"process_range memory is fixed ({memory[0]!r}): retrying an OOM at the "
        "same request only burns the retries"
    )


def test_range_memory_is_a_memory_literal_everywhere():
    string_valued = []

    for config in CONFIGS:
        if ".venv" in config.parts or "work" in config.parts:
            continue
        for value in RANGE_MEMORY.findall(config.read_text()):
            if value.startswith(("'", '"')):
                string_valued.append(f"{config.relative_to(ROOT)}: {value}")

    assert string_valued == [], (
        "range.memory is multiplied by task.attempt, so a string value is "
        f"repeated rather than scaled: {string_valued}"
    )
