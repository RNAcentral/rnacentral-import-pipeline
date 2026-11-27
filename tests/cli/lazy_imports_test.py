# -*- coding: utf-8 -*-

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_importing_the_cli_does_not_import_matplotlib():
    # Every rnac task pays for matplotlib's font cache otherwise, and under --contain it warns.
    code = "import sys, rnacentral_pipeline.cli; sys.exit('matplotlib' in sys.modules)"
    result = subprocess.run([sys.executable, "-c", code], cwd=ROOT)

    assert result.returncode == 0, (
        "importing rnacentral_pipeline.cli pulls in matplotlib"
    )
