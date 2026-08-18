import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARSE_DATABASES_WORKFLOW = ROOT / "workflows" / "parse-databases.nf"

DATABASE_INCLUDE_RE = re.compile(
    r"^include\s+\{[^}]+\}\s+from\s+'\.\/databases\/([^']+)'",
    re.MULTILINE,
)
# A workflow is guarded either by a process level `when:` or by wrapping the
# whole body in `if ( params.databases.<name>?.run )`. Both keep the work from
# running; the second skips the process entirely.
PROCESS_RUN_GUARD_RE = re.compile(
    r"^\s*when:\s*(?:\{[^\n]*params\.databases[^\n]*\.run[^\n]*\}|"
    r"params\.databases[^\n]*\.run)\s*$"
    r"|^\s*when:\s*\n\s*params\.databases[^\n]*\.run",
    re.MULTILINE,
)
WORKFLOW_RUN_GUARD_RE = re.compile(
    r"^\s*if\s*\(\s*params\.databases[^\n]*\.run\s*\)",
    re.MULTILINE,
)


def is_guarded(text):
    return bool(PROCESS_RUN_GUARD_RE.search(text) or WORKFLOW_RUN_GUARD_RE.search(text))


UNSAFE_DATABASE_RUN_RE = re.compile(
    r"params\.databases(?:\.[A-Za-z_][A-Za-z0-9_]*|\[[^\]]+\])+" r"(?<!\?)\.run"
)


def database_workflows():
    modules = DATABASE_INCLUDE_RE.findall(PARSE_DATABASES_WORKFLOW.read_text())
    assert modules
    return [ROOT / "workflows" / "databases" / f"{module}.nf" for module in modules]


def test_included_database_workflows_have_run_guards():
    missing = [
        str(path.relative_to(ROOT))
        for path in database_workflows()
        if not is_guarded(path.read_text())
    ]

    assert missing == []


def test_database_run_checks_are_null_safe():
    unsafe = []

    for path in database_workflows():
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            if UNSAFE_DATABASE_RUN_RE.search(line):
                unsafe.append(f"{path.relative_to(ROOT)}:{lineno}")

    assert unsafe == []


def test_taxonomy_context_uses_global_config_flag():
    workflow = PARSE_DATABASES_WORKFLOW.read_text()

    assert "if (params.get('needs_taxonomy', false))" in workflow
    assert "build_context | set { context }" in workflow
    assert "channel.empty() | set { context }" in workflow
