import importlib.util
import io
from pathlib import Path

import numpy as np
import polars as pl

BIN = Path(__file__).resolve().parents[2] / "bin" / "litscan-scan-job.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("litscan_scan_job", BIN)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


scan_job = _load_script()


class FakePipeline:
    def predict_proba(self, texts):
        return np.array([[0.7, 0.3]])


def test_rna_related_writes_as_boolean():
    # A numpy bool in a row dict became Float64 in polars, writing "0.0"
    # into a boolean column and failing the load.
    probability, rna_related = scan_job.classify_abstract("text", FakePipeline())
    buf = io.StringIO()
    pl.DataFrame([{"rna_related": rna_related}]).write_csv(buf, include_header=False)
    assert buf.getvalue().strip() == "false"
