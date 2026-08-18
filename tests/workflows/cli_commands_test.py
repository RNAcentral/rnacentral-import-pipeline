# -*- coding: utf-8 -*-

"""
A workflow reaches its parser through one string: the ``rnac <command>`` in its
script block. Nothing at import time connects the two, so a parser can sit in
the tree for years with no way to invoke it -- which is exactly what happened to
MGI, ported from luigi in 2017 and never given a CLI command or a workflow.

These check both halves of that wiring: the command a workflow calls exists, and
the module defining it is actually registered on the top level CLI.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLI_INIT = ROOT / "rnacentral_pipeline" / "cli" / "__init__.py"

# Groups live either in the cli package or beside the database they belong to.
CLI_FILES = sorted(
    list((ROOT / "rnacentral_pipeline" / "cli").glob("*.py"))
    + list((ROOT / "rnacentral_pipeline" / "databases").glob("*/cli.py"))
)

DEFINITION_RE = re.compile(r'@click\.(?:group|command)\(\s*"([^"]+)"')
INVOCATION_RE = re.compile(r"\brnac\s+([a-z0-9][a-z0-9_-]*)")
REGISTRATION_RE = re.compile(r"^cli\.add_command\((\w+)\.", re.MULTILINE)


def definitions():
    """Map every CLI command name to the module token that must be registered."""
    found = {}
    for path in CLI_FILES:
        # databases/<name>/cli.py is imported as <name>, cli/<name>.py as <name>.
        token = path.parent.name if path.stem == "cli" else path.stem
        for name in DEFINITION_RE.findall(path.read_text()):
            found[name] = token
    return found


def invocations():
    """Map every rnac command a workflow calls to the files calling it."""
    found = {}
    for path in sorted(ROOT.glob("workflows/**/*.nf")):
        for name in INVOCATION_RE.findall(path.read_text()):
            found.setdefault(name, set()).add(str(path.relative_to(ROOT)))
    assert found
    return found


def test_every_command_a_workflow_calls_exists():
    known = definitions()
    missing = {
        name: sorted(paths)
        for name, paths in invocations().items()
        if name not in known
    }
    assert missing == {}


def test_every_command_a_workflow_calls_is_registered():
    known = definitions()
    registered = set(REGISTRATION_RE.findall(CLI_INIT.read_text()))
    assert registered

    unregistered = {
        name: known[name]
        for name in invocations()
        if name in known and known[name] not in registered
    }
    assert unregistered == {}
