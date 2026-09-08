#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
One-off conversion: gensim Word2Vec .emb -> plain parquet embedding table
(so_term: str, vector: list[float32]).

gensim is not a pipeline dependency any more (no Python 3.14 wheel, and the
pipeline only ever did KeyedVectors lookups on this file). Rerun this
whenever the SO Node2Vec model is retrained, with gensim installed just for
the conversion:

    uv run --with gensim python bin/convert-so-model-to-parquet.py \
        so_embedding_model.emb so_embedding_model.parquet
"""

import argparse

import polars as pl
from gensim.models import Word2Vec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model_path", help="gensim Word2Vec .emb model file")
    parser.add_argument("output_parquet")
    args = parser.parse_args()

    model = Word2Vec.load(args.model_path)
    terms = list(model.wv.key_to_index)
    vectors = [model.wv[term].tolist() for term in terms]

    pl.DataFrame({"so_term": terms, "vector": vectors}).write_parquet(
        args.output_parquet
    )
    print(
        f"Wrote {len(terms)} vectors ({model.wv.vector_size}-d) to {args.output_parquet}"
    )


if __name__ == "__main__":
    main()
