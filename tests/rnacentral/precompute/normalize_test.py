# -*- coding: utf-8 -*-

import polars as pl

from rnacentral_pipeline.rnacentral.precompute import normalize


def _accession(**overrides):
    base = {
        "id": 1,
        "urs_id": 1,
        "urs_taxid": "URS1_9606",
        "accession": "A1",
        "last_release": 5,
        "is_active": True,
        "description": "d",
        "gene": None,
        "optional_id": None,
        "database": "ENA",
        "species": None,
        "common_name": None,
        "feature_name": None,
        "ncrna_class": None,
        "locus_tag": None,
        "organelle": None,
        "lineage": None,
        "all_species": ["Homo sapiens", None],
        "all_common_names": [None, None],
        "so_rna_type": "SO:0000252",
    }
    base.update(overrides)
    return base


def _write_accessions_parquet(path, rows):
    pl.DataFrame(rows).write_parquet(path)


def test_read_accessions_groups_raw_rows_by_id(tmp_path):
    # get-accessions/query.sql's actual output shape: one flat row per
    # accession, ordered by id but with no grouping - an id can have
    # multiple rows (multiple xrefs to the same urs_taxid).
    path = tmp_path / "accessions.parquet"
    _write_accessions_parquet(
        path,
        [
            _accession(id=1, accession="A1"),
            _accession(id=1, accession="A1b"),
            _accession(id=2, accession="A2"),
        ],
    )

    df = normalize.read_accessions(path)

    by_id = {row["id"]: row["data"] for row in df.iter_rows(named=True)}
    assert set(by_id.keys()) == {1, 2}
    assert {a["accession"] for a in by_id[1]} == {"A1", "A1b"}
    assert [a["accession"] for a in by_id[2]] == ["A2"]


def _metadata_row(id_: int = 1, **overrides) -> dict:
    row = {
        "id": id_,
        "urs_id": id_,
        "urs_taxid": f"URS{id_}_9606",
        "urs": f"URS{id_}",
        "taxid": 9606,
        "length": 100,
        "coordinates": [],
        "previous": None,
        "rfam_hits": [],
        "r2dt_hits": None,
        "orf_info": None,
        "possible_orf": None,
        "possible_orf_stopfree": None,
        "possible_orf_tcode": None,
    }
    row.update(overrides)
    return row


def test_read_metadata_reads_nested_structs(tmp_path):
    path = tmp_path / "metadata.parquet"
    pl.DataFrame(
        [_metadata_row(1), _metadata_row(2, orf_info={"sources": ["cpat"]})]
    ).write_parquet(path)

    df = normalize.read_metadata(path)

    assert df["urs"].to_list() == ["URS1", "URS2"]
    assert df["orf_info"].to_list() == [None, {"sources": ["cpat"]}]


def test_join_accessions_metadata_is_an_inner_join_on_id():
    accessions_df = pl.DataFrame({"id": [1, 2], "data": [["a"], ["b"]]})
    metadata_df = pl.DataFrame({"id": [2, 3], "urs": ["URS2", "URS3"]})

    joined = normalize.join_accessions_metadata(accessions_df, metadata_df)

    # id=1 (accessions only) and id=3 (metadata only) are both dropped -
    # only id=2, present on both sides, survives.
    assert joined["id"].to_list() == [2]
    assert joined["urs"].to_list() == ["URS2"]


def _joined_row(**overrides):
    base = {
        "data": [_accession()],
        "urs": "URS1",
        "taxid": 9606,
        "length": 100,
        "coordinates": [],
        "previous": None,
        "rfam_hits": [],
        "r2dt_hits": None,
        "orf_info": None,
        "possible_orf": None,
        "possible_orf_stopfree": None,
        "possible_orf_tcode": None,
    }
    base.update(overrides)
    return base


def _normalize(row: dict) -> dict:
    return normalize.normalize(pl.DataFrame([row])).row(0, named=True)


def test_normalize_row_drops_null_entries_from_species_and_common_names():
    row = _joined_row()

    result = _normalize(row)

    assert result["accessions"][0]["all_species"] == ["Homo sapiens"]
    assert result["accessions"][0]["all_common_names"] == []


def test_normalize_row_computes_last_release_and_deleted_over_whole_group():
    row = _joined_row(
        data=[
            _accession(last_release=5, is_active=False),
            _accession(last_release=9, is_active=True),
            _accession(last_release=3, is_active=False),
        ]
    )

    result = _normalize(row)

    assert result["last_release"] == 9
    assert result["deleted"] is False  # one accession is still active


def test_normalize_row_deleted_true_when_all_accessions_inactive():
    row = _joined_row(data=[_accession(is_active=False), _accession(is_active=False)])

    result = _normalize(row)

    assert result["deleted"] is True


def test_normalize_row_wraps_present_r2dt_hit_as_single_element_list():
    row = _joined_row(r2dt_hits={"model_id": 5540, "model_name": "EC_SSU_3D"})

    result = _normalize(row)

    assert result["r2dt_hits"] == [{"model_id": 5540, "model_name": "EC_SSU_3D"}]


def test_normalize_row_r2dt_hits_is_empty_list_when_absent():
    row = _joined_row(r2dt_hits=None)

    result = _normalize(row)

    assert result["r2dt_hits"] == []


def test_normalize_row_output_has_exactly_the_normalized_fields():
    row = _joined_row()

    result = _normalize(row)

    assert set(result.keys()) == {
        "urs",
        "taxid",
        "length",
        "last_release",
        "coordinates",
        "accessions",
        "deleted",
        "previous",
        "rfam_hits",
        "r2dt_hits",
        "orf_info",
        "possible_orf",
        "possible_orf_stopfree",
        "possible_orf_tcode",
    }


def test_write_reads_joins_normalizes_and_writes(tmp_path):
    accessions_path = tmp_path / "accessions.parquet"
    metadata_path = tmp_path / "metadata.parquet"
    output_path = tmp_path / "output.parquet"

    _write_accessions_parquet(
        accessions_path,
        [
            _accession(id=1),
            # id=2 has no accession rows at all - never appears as a group.
            _accession(id=3),  # no metadata match
        ],
    )
    pl.DataFrame([_metadata_row(1)]).write_parquet(metadata_path)

    normalize.write(accessions_path, metadata_path, output_path)

    out = pl.read_parquet(output_path)
    assert out.height == 1
    result = out.row(0, named=True)
    assert result["urs"] == "URS1"
