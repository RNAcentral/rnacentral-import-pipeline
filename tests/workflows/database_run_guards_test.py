from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
PARSE_DATABASES_WORKFLOW = ROOT / "workflows" / "parse-databases.nf"

DATABASE_INCLUDE_RE = re.compile(
    r"^include\s+\{[^}]+\}\s+from\s+'\.\/databases\/([^']+)'",
    re.MULTILINE,
)
DATABASE_RUN_GUARD_RE = re.compile(
    r"^\s*when:\s*(?:\{[^\n]*params\.databases[^\n]*\.run[^\n]*\}|"
    r"params\.databases[^\n]*\.run)\s*$"
    r"|^\s*when:\s*\n\s*params\.databases[^\n]*\.run",
    re.MULTILINE,
)


def database_workflows():
    modules = DATABASE_INCLUDE_RE.findall(PARSE_DATABASES_WORKFLOW.read_text())
    assert modules
    return [ROOT / "workflows" / "databases" / f"{module}.nf" for module in modules]


def test_included_database_workflows_have_run_guards():
    missing = [
        str(path.relative_to(ROOT))
        for path in database_workflows()
        if not DATABASE_RUN_GUARD_RE.search(path.read_text())
    ]

    assert missing == []
