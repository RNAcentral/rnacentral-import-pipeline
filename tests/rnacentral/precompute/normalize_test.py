# -*- coding: utf-8 -*-

import json

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


def test_read_metadata_reads_nested_structs(tmp_path):
    path = tmp_path / "metadata.json"
    path.write_text(
        '{"id":1,"urs_id":1,"urs_taxid":"URS1_9606","urs":"URS1","taxid":9606,'
        '"length":100,"coordinates":[],"previous":null,"rfam_hits":[],'
        '"r2dt_hits":null,"orf_info":null,"possible_orf":null,'
        '"possible_orf_stopfree":null,"possible_orf_tcode":null}\n'
    )

    df = normalize.read_metadata(path)

    assert df["urs"].to_list() == ["URS1"]
    assert df["r2dt_hits"].to_list() == [None]


def _metadata_line(id_: int, **overrides) -> str:
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
    return json.dumps(row)


def test_read_metadata_handles_a_column_that_is_null_past_the_inference_sample(
    tmp_path,
):
    # pl.read_ndjson infers each column's type from only the first 100 rows
    # by default. orf_info is null for the vast majority of real rows and
    # only non-null occasionally (CPAT-flagged sequences) - a real
    # metadata.json can easily have hundreds of thousands of null orf_info
    # rows before the first real one. Reproduce that shape at a size a test
    # can run fast: 150 null rows (past the 100-row default sample), then
    # one row with a real orf_info value.
    path = tmp_path / "metadata.json"
    lines = [_metadata_line(i) for i in range(1, 151)]
    lines.append(_metadata_line(151, orf_info={"sources": ["cpat"]}))
    path.write_text("\n".join(lines) + "\n")

    df = normalize.read_metadata(path)

    assert df.height == 151
    assert df.filter(pl.col("id") == 151)["orf_info"].to_list() == [
        {"sources": ["cpat"]}
    ]


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


def test_normalize_row_drops_null_entries_from_species_and_common_names():
    row = _joined_row()

    result = normalize.normalize_row(row)

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

    result = normalize.normalize_row(row)

    assert result["last_release"] == 9
    assert result["deleted"] is False  # one accession is still active


def test_normalize_row_deleted_true_when_all_accessions_inactive():
    row = _joined_row(data=[_accession(is_active=False), _accession(is_active=False)])

    result = normalize.normalize_row(row)

    assert result["deleted"] is True


def test_normalize_row_wraps_present_r2dt_hit_as_single_element_list():
    row = _joined_row(r2dt_hits={"model_id": 5540, "model_name": "EC_SSU_3D"})

    result = normalize.normalize_row(row)

    assert result["r2dt_hits"] == [{"model_id": 5540, "model_name": "EC_SSU_3D"}]


def test_normalize_row_r2dt_hits_is_empty_list_when_absent():
    row = _joined_row(r2dt_hits=None)

    result = normalize.normalize_row(row)

    assert result["r2dt_hits"] == []


def test_normalize_row_output_has_exactly_the_normalized_fields():
    row = _joined_row()

    result = normalize.normalize_row(row)

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


def test_write_output_writes_one_json_object_per_line(tmp_path):
    path = tmp_path / "out.json"

    count = normalize.write_output([{"a": 1}, {"a": 2}], path)

    assert count == 2
    lines = path.read_text().splitlines()
    assert [json.loads(line) for line in lines] == [{"a": 1}, {"a": 2}]


def test_write_reads_joins_normalizes_and_writes(tmp_path):
    accessions_path = tmp_path / "accessions.parquet"
    metadata_path = tmp_path / "metadata.json"
    output_path = tmp_path / "output.json"

    _write_accessions_parquet(
        accessions_path,
        [
            _accession(id=1),
            # id=2 has no accession rows at all - never appears as a group.
            _accession(id=3),  # no metadata match
        ],
    )
    metadata_row = {k: v for k, v in _joined_row().items() if k != "data"}
    metadata_row["id"] = 1
    metadata_row["urs_id"] = 1
    metadata_row["urs_taxid"] = "URS1_9606"
    metadata_path.write_text(json.dumps(metadata_row) + "\n")

    normalize.write(accessions_path, metadata_path, output_path)

    lines = output_path.read_text().splitlines()
    assert len(lines) == 1
    result = json.loads(lines[0])
    assert result["urs"] == "URS1"
