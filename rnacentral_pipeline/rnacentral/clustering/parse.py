# -*- coding: utf-8 -*-

"""
Flatten the parquet output of rnacentral-mmseqs-clustering into the shapes of
load_sequence_cluster and load_sequence_absorption (files/clustering/schema.sql).
"""

from pathlib import Path

import polars as pl

# Which member file holds which kind. The producer splits on its
# fragment_min_length_ratio, so the file is the only record of the kind.
ABSORBED_FILES = {"fragments": "fragment", "variants": "variant"}


def parse(src: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)

    # clusters.parquet also carries length, so_types and the is_representative_* flags;
    # all derivable from rna / rnc_rna_precomputed / the ids, so they are not loaded.
    pl.scan_parquet(src / "clusters.parquet").select(
        "urs", "cluster_95", "cluster_90"
    ).sink_parquet(out / "clusters.parquet")

    frames = [
        pl.scan_parquet(path)
        .select("parent_urs", "fragment_urs", "fragment_identity")
        .explode("fragment_urs", "fragment_identity")
        .select(
            urs="fragment_urs",
            parent_urs="parent_urs",
            absorbed_as=pl.lit(kind),
            identity=pl.col("fragment_identity").cast(pl.Float32),
        )
        for name, kind in ABSORBED_FILES.items()
        if (path := src / f"{name}.parquet").exists()
    ]
    if not frames:
        raise ValueError(f"neither fragments.parquet nor variants.parquet in {src}")
    pl.concat(frames).sink_parquet(out / "absorption.parquet")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        src, out = Path(tmp, "src"), Path(tmp, "out")
        src.mkdir()
        pl.DataFrame(
            {
                "urs": ["A", "B", "C", "D"],
                "cluster_95": ["A", "A", "A", "D"],
                "cluster_90": ["A", "A", "A", "A"],
            }
        ).write_parquet(src / "clusters.parquet")
        member = {
            "parent_urs": ["A"],
            "parent_len": [100],
            "fragment_len": [[50]],
            "fragment_ratio": [[0.5]],
            "n_fragments": [1],
        }
        pl.DataFrame(
            {**member, "fragment_urs": [["B"]], "fragment_identity": [[0.97]]}
        ).write_parquet(src / "fragments.parquet")
        pl.DataFrame(
            {**member, "fragment_urs": [["C"]], "fragment_identity": [[None]]},
            schema_overrides={"fragment_identity": pl.List(pl.Float32)},
        ).write_parquet(src / "variants.parquet")

        parse(src, out)
        absorbed = pl.read_parquet(out / "absorption.parquet").sort("urs")
        assert absorbed["urs"].to_list() == ["B", "C"], absorbed
        assert absorbed["absorbed_as"].to_list() == ["fragment", "variant"], absorbed
        assert absorbed["identity"].null_count() == 1, absorbed
        assert pl.read_parquet(out / "clusters.parquet").columns == [
            "urs",
            "cluster_95",
            "cluster_90",
        ]
    print("ok")
