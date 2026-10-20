# -*- coding: utf-8 -*-

"""
Tests for the parquet-based SO embedding loader, which replaced loading a
gensim Word2Vec .emb file directly (gensim was dropped -- no Python 3.14
wheel, and the pipeline only ever did KeyedVectors lookups on it).
"""

import sys
import types

import numpy as np
import polars as pl

# preprocessing.py imports the Rust gene_preprocessing extension at module
# level for unrelated functions; load_so_model doesn't touch it, so stub it
# out just long enough to import -- and remove the stub right after, so we
# don't fool other test modules' `pytest.importorskip("gene_preprocessing")`
# into thinking the real extension is installed.
_stubbed_gene_preprocessing = False
if "gene_preprocessing" not in sys.modules:
    try:
        import gene_preprocessing  # noqa: F401
    except ImportError:
        sys.modules["gene_preprocessing"] = types.ModuleType("gene_preprocessing")
        _stubbed_gene_preprocessing = True

from rnacentral_pipeline.rnacentral.genes.random_forest.preprocessing import (
    load_so_model,
)

if _stubbed_gene_preprocessing:
    del sys.modules["gene_preprocessing"]


def test_load_so_model_normalises_vectors(tmp_path):
    path = tmp_path / "so_model.parquet"
    pl.DataFrame(
        {
            "so_term": ["SO:A", "SO:B"],
            "vector": [[3.0, 4.0], [1.0, 0.0]],
        }
    ).write_parquet(path)

    model = load_so_model(str(path))

    assert set(model) == {"SO:A", "SO:B"}
    np.testing.assert_allclose(model["SO:A"], [0.6, 0.8], atol=1e-6)
    np.testing.assert_allclose(float(np.linalg.norm(model["SO:B"])), 1.0)
