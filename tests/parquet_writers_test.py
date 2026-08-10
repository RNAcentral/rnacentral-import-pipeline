# -*- coding: utf-8 -*-

"""
Unit tests for the shared streaming Parquet writer helpers.

Two things are pinned down here:

1. The row-length guards on ``ParquetTable`` / ``TypedParquetWrapper``. Without
   them a short row silently truncates and a long row silently drops columns.

2. The str -> typed conversion contract. Every parser in this repo still emits
   stringly-typed rows from its CSV-era ``writeable()`` methods, while the
   schemas in :mod:`rnacentral_pipeline.schemas` are typed. Getting that bridge
   wrong is what produced the class of runtime crashes this migration hit --
   e.g. ``Reference.writeable()`` yielding ``""`` for a missing pmid against an
   ``int64`` column.
"""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from rnacentral_pipeline import schemas
from rnacentral_pipeline.parquet_writers import (
    TypedParquetWrapper,
    converter_for,
    parquet_writer,
    typed_parquet_writer,
)

SIMPLE = pa.schema(
    [
        pa.field("name", pa.string(), nullable=False),
        pa.field("count", pa.int64()),
    ]
)


# ---------------------------------------------------------------------------
# Row length guards


def test_writerow_rejects_short_row(tmp_path):
    with parquet_writer(tmp_path / "out.parquet", SIMPLE) as table:
        with pytest.raises(ValueError, match="Row length"):
            table.writerow(("only-one",))


def test_writerow_rejects_long_row(tmp_path):
    with parquet_writer(tmp_path / "out.parquet", SIMPLE) as table:
        with pytest.raises(ValueError, match="Row length"):
            table.writerow(("a", 1, "extra"))


def test_writerows_rejects_bad_row(tmp_path):
    with parquet_writer(tmp_path / "out.parquet", SIMPLE) as table:
        with pytest.raises(ValueError, match="Row length"):
            table.writerows([("a", 1), ("b",)])


def test_typed_wrapper_rejects_long_row(tmp_path):
    """A long row must raise rather than have zip() silently drop the tail."""
    with parquet_writer(tmp_path / "out.parquet", SIMPLE) as table:
        wrapper = TypedParquetWrapper(table, SIMPLE)
        with pytest.raises(ValueError, match="Row length"):
            wrapper.writerow(("a", "1", "extra"))


# ---------------------------------------------------------------------------
# Converters


@pytest.mark.parametrize("blank", ["", "NaN", "nan", None])
def test_nullable_int_converter_maps_blanks_to_none(blank):
    convert = converter_for(pa.field("x", pa.int64(), nullable=True))
    assert convert(blank) is None


@pytest.mark.parametrize("blank", ["", "NaN", "nan", None])
def test_nullable_float_converter_maps_blanks_to_none(blank):
    convert = converter_for(pa.field("x", pa.float64(), nullable=True))
    assert convert(blank) is None


def test_non_nullable_int_converter_raises_on_blank():
    convert = converter_for(pa.field("x", pa.int64(), nullable=False))
    with pytest.raises(ValueError):
        convert("")


@pytest.mark.parametrize(
    "value,expected",
    [
        ("1", True),
        ("t", True),
        ("true", True),
        ("True", True),
        (True, True),
        ("0", False),
        ("f", False),
        ("false", False),
        ("False", False),
        (False, False),
    ],
)
def test_bool_converter_accepts_csv_and_native_tokens(value, expected):
    convert = converter_for(pa.field("x", pa.bool_(), nullable=False))
    assert convert(value) is expected


def test_non_nullable_bool_converter_raises_on_nonsense():
    convert = converter_for(pa.field("x", pa.bool_(), nullable=False))
    with pytest.raises(ValueError):
        convert("maybe")


def test_nullable_bool_converter_maps_blank_to_none():
    convert = converter_for(pa.field("x", pa.bool_(), nullable=True))
    assert convert("") is None


def test_string_converter_is_identity():
    convert = converter_for(pa.field("x", pa.string()))
    assert convert("anything") == "anything"
    assert convert(None) is None


# ---------------------------------------------------------------------------
# Round trips


def test_typed_writer_round_trip(tmp_path):
    out = tmp_path / "out.parquet"
    with typed_parquet_writer(out, SIMPLE) as writer:
        writer.writerows([("a", "1"), ("b", "")])

    table = pq.read_table(out)
    assert table.schema.equals(SIMPLE)
    assert table.to_pylist() == [
        {"name": "a", "count": 1},
        {"name": "b", "count": None},
    ]


def test_writer_flushes_across_batches(tmp_path):
    """Rows must survive a mid-stream flush, not just the final close()."""
    out = tmp_path / "out.parquet"
    with parquet_writer(out, SIMPLE, batch_size=2) as table:
        table.writerows([("r%d" % i, i) for i in range(5)])

    result = pq.read_table(out)
    assert result.num_rows == 5
    assert result.column("name").to_pylist() == ["r0", "r1", "r2", "r3", "r4"]


def test_writer_creates_parent_directories(tmp_path):
    out = tmp_path / "nested" / "deeper" / "out.parquet"
    with parquet_writer(out, SIMPLE) as table:
        table.writerow(("a", 1))
    assert out.is_file()


# ---------------------------------------------------------------------------
# Schema conformance across every declared schema.
#
# This is the broad guard for the migration: it drives a synthetic stringly
# typed row -- the shape every CSV-era ``writeable()`` still emits -- through
# TypedParquetWrapper for each schema and asserts the file that lands matches
# the declared schema exactly. A future edit that reintroduces a stringified
# int, or a schema whose types drift from what the writers emit, fails here
# rather than during a production load.

ALL_SCHEMAS = sorted(
    (name, schema)
    for name, schema in vars(schemas).items()
    if isinstance(schema, pa.Schema)
)


def _stringly_value(field: pa.Field, blank: bool):
    """
    A value in the form a CSV-era ``writeable()`` would emit for ``field``.

    When ``blank`` and the field is nullable we return ``""`` -- the empty
    field a csv.writer produces for a missing value, and the exact input that
    broke the untyped writers.
    """
    if blank and field.nullable:
        return ""
    t = field.type
    if pa.types.is_boolean(t):
        return "1"
    if pa.types.is_integer(t):
        return "42"
    if pa.types.is_floating(t):
        return "1.5"
    if pa.types.is_list(t):
        # Array columns arrive as the pg array literal writeable() emits.
        return '{"a","b"}'
    return 'some,text with "quotes"'


@pytest.mark.parametrize("name,schema", ALL_SCHEMAS, ids=[n for n, _ in ALL_SCHEMAS])
def test_schema_accepts_populated_stringly_row(tmp_path, name, schema):
    out = tmp_path / f"{name}.parquet"
    with typed_parquet_writer(out, schema) as writer:
        writer.writerow([_stringly_value(f, blank=False) for f in schema])

    table = pq.read_table(out)
    assert table.schema.equals(schema)
    assert table.num_rows == 1


@pytest.mark.parametrize("name,schema", ALL_SCHEMAS, ids=[n for n, _ in ALL_SCHEMAS])
def test_schema_accepts_blank_nullables(tmp_path, name, schema):
    """
    Every nullable column must accept the empty string a csv.writer emits for a
    missing value, and land as NULL rather than crashing or storing "".
    """
    out = tmp_path / f"{name}.parquet"
    with typed_parquet_writer(out, schema) as writer:
        writer.writerow([_stringly_value(f, blank=True) for f in schema])

    table = pq.read_table(out)
    assert table.schema.equals(schema)

    row = table.to_pylist()[0]
    for field in schema:
        if not field.nullable:
            continue
        value = row[field.name]
        if pa.types.is_string(field.type):
            assert value == ""
        else:
            assert value is None, f"{name}.{field.name} kept {value!r}"


@pytest.mark.parametrize("name,schema", ALL_SCHEMAS, ids=[n for n, _ in ALL_SCHEMAS])
def test_schema_rejects_wrong_width_row(tmp_path, name, schema):
    out = tmp_path / f"{name}.parquet"
    with typed_parquet_writer(out, schema) as writer:
        row = [_stringly_value(f, blank=False) for f in schema]
        with pytest.raises(ValueError, match="Row length"):
            writer.writerow(row + ["extra"])
        # Leave a valid row so the file is not empty on close.
        writer.writerow(row)


def test_no_schema_is_entirely_non_nullable_by_accident():
    """
    Sanity check on the blank-nullable test above: it only has teeth if some
    schemas actually declare nullable columns. If a refactor ever marked
    everything NOT NULL this would go green while testing nothing.
    """
    nullable_counts = {
        name: sum(1 for f in schema if f.nullable) for name, schema in ALL_SCHEMAS
    }
    assert sum(nullable_counts.values()) > 0
