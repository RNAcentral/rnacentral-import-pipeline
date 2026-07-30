# -*- coding: utf-8 -*-

"""
End-to-end conformance for the parquet writers, driven by real domain objects.

``tests/parquet_writers_test.py`` checks the schemas accept synthetic rows.
This checks the other half: that the values our actual ``writeable()`` methods
produce are accepted by those schemas. Every one of the type-mismatch crashes
this migration hit lived in that gap -- a producer emitting a stringified int
against an int64 column -- and it only surfaced at load time.

The awkward cases are deliberate: a reference with no pmid (the empty string
that broke ``REFERENCES``), a description containing a comma and a quote, and
an entry with none of the optional collections populated.
"""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from rnacentral_pipeline import schemas, writers
from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.data import regions


def _entry(**kwargs):
    defaults = dict(
        primary_id="primary",
        accession="ACC.1",
        ncbi_tax_id=9606,
        database="ENA",
        sequence="ACGTACGTACGTACGT",
        regions=[],
        rna_type="SO:0000655",
        url="https://example.org/ACC.1",
        seq_version="1",
        species="Homo sapiens",
        common_name="human",
        lineage="Eukaryota; Metazoa; Homo sapiens",
        description='A ncRNA, with a comma and a " quote',
        gene="GENE1",
        product="a product",
        mol_type="genomic DNA",
    )
    defaults.update(kwargs)
    return data.Entry(**defaults)


def _reference(pmid=12345, doi="10.1000/xyz"):
    return data.Reference(
        authors="Smith, J. and Jones, A.",
        location="Journal of Things, 1(2)",
        title='A study of "things"',
        pmid=pmid,
        doi=doi,
    )


def _region():
    return regions.SequenceRegion(
        assembly_id="GRCh38",
        chromosome="1",
        strand=1,
        exons=[regions.Exon(start=100, stop=200)],
        coordinate_system=regions.CoordinateSystem.zero_based(),
    )


def _go_annotation():
    return data.GoTermAnnotation(
        rna_id="URS0000000001_9606",
        qualifier="part_of",
        term_id="GO:0003723",
        evidence_code="ECO:0000305",
        extensions=[],
        assigned_by="RNAcentral",
        publications=[data.IdReference.build(12345)],
    )


def _written_tables(path):
    """Every parquet file the writer produced, keyed by logical name."""
    return {p.stem: pq.read_table(p) for p in sorted(path.glob("*.parquet"))}


ENTRY_CASES = {
    "minimal": lambda: _entry(),
    "with_references": lambda: _entry(references=[_reference()]),
    "reference_without_pmid": lambda: _entry(references=[_reference(pmid=None)]),
    "reference_without_doi": lambda: _entry(references=[_reference(doi=None)]),
    "with_regions": lambda: _entry(regions=[_region()]),
    "with_go_annotations": lambda: _entry(go_annotations=[_go_annotation()]),
    "long_sequence": lambda: _entry(sequence="ACGT" * 1000),
    "no_optional_metadata": lambda: _entry(
        species=None,
        common_name=None,
        lineage=None,
        gene=None,
        product=None,
        mol_type=None,
    ),
    "with_secondary_structure": lambda: _entry(
        secondary_structure=data.SecondaryStructure(dot_bracket="((((....))))")
    ),
}

# Pre-existing bug, not introduced by this PR and present on dev too:
# SecondaryStructure.md5 passes a str to helpers.hashes.md5, which calls
# hashlib.md5() and needs bytes. Any entry carrying a secondary structure
# raises TypeError, on the CSV path as much as the parquet one. xfail_strict
# is on repo-wide, so this turns red the moment someone fixes it.
BROKEN_UPSTREAM = {
    "with_secondary_structure": pytest.mark.xfail(
        reason="hashes.md5 needs bytes; SecondaryStructure.md5 passes str",
        raises=TypeError,
        strict=True,
    )
}

ENTRY_PARAMS = [
    pytest.param(name, marks=[BROKEN_UPSTREAM[name]] if name in BROKEN_UPSTREAM else [])
    for name in sorted(ENTRY_CASES)
]


@pytest.mark.parametrize("case", ENTRY_PARAMS)
def test_entry_writer_output_matches_declared_schemas(tmp_path, case):
    entry = ENTRY_CASES[case]()

    with writers.parquet_entry_writer(tmp_path) as writer:
        writer.write([entry])

    written = _written_tables(tmp_path)
    assert written, "writer produced no parquet files"

    for name, table in written.items():
        expected = schemas.ENTRY_WRITER_SCHEMAS[name]
        assert table.schema.equals(expected), f"{name} schema drifted"


def test_entry_writer_drops_tables_that_got_no_rows(tmp_path):
    """
    A table with no rows must leave no file behind.

    Nothing reads these by a fixed name -- every consumer globs, and
    load-data.nf builds its channel from the files that exist -- so an empty
    parquet file only buys a merge_and_import task that loads nothing. The CSV
    path gets this for free: an unused table is a 0-byte file, which
    load-data.nf's `filter { f -> !f.isEmpty() }` already drops. A zero-row
    parquet file is a few hundred bytes of header and footer, so it does not.
    """
    with writers.parquet_entry_writer(tmp_path) as writer:
        writer.write([_entry()])

    written = _written_tables(tmp_path)
    assert written, "writer produced no parquet files"
    assert set(written) <= set(schemas.ENTRY_WRITER_SCHEMAS)
    for name, table in written.items():
        assert table.num_rows, f"{name}.parquet was written with no rows"

    # The entry fixture populates some tables and not others; if that ever
    # stops being true this test is no longer testing anything.
    assert set(written) != set(schemas.ENTRY_WRITER_SCHEMAS)


def test_entry_writer_round_trips_awkward_text(tmp_path):
    """Commas and quotes must survive as data, not become structure."""
    description = 'A ncRNA, with a comma and a " quote'
    with writers.parquet_entry_writer(tmp_path) as writer:
        writer.write([_entry(description=description)])

    accessions = pq.read_table(tmp_path / "accessions.parquet").to_pylist()
    assert any(description in str(v) for row in accessions for v in row.values())


def test_reference_without_pmid_lands_as_null_not_empty_string(tmp_path):
    """
    REFERENCES declares pmid as int64 while Reference.writeable() stringifies.
    A missing pmid must become NULL rather than crash on "" -- this is the
    exact failure the typed writers were introduced to prevent.
    """
    with writers.parquet_entry_writer(tmp_path) as writer:
        writer.write([_entry(references=[_reference(pmid=None)])])

    table = pq.read_table(tmp_path / "references.parquet")
    assert table.schema.equals(schemas.REFERENCES)
    if table.num_rows:
        pmids = table.column("pmid").to_pylist()
        assert all(p is None or isinstance(p, int) for p in pmids)


def test_writer_raises_when_given_no_entries(tmp_path):
    with pytest.raises(ValueError):
        with writers.parquet_entry_writer(tmp_path) as writer:
            writer.write([])


def test_writer_closes_files_when_write_raises(tmp_path):
    """
    An exception mid-write must still leave readable parquet files rather than
    truncated ones holding an open handle.
    """
    with pytest.raises(ValueError):
        with writers.parquet_entry_writer(tmp_path) as writer:
            writer.write([])

    for path in tmp_path.glob("*.parquet"):
        pq.read_schema(path)


def test_multiple_entries_accumulate(tmp_path):
    entries = [_entry(accession=f"ACC.{i}", primary_id=f"p{i}") for i in range(10)]
    with writers.parquet_entry_writer(tmp_path) as writer:
        writer.write(entries)

    accessions = pq.read_table(tmp_path / "accessions.parquet")
    assert accessions.num_rows == len(entries)


def test_entry_writer_schemas_are_all_arrow_schemas():
    """Guards against a plain dict or list sneaking into the registry."""
    for name, schema in schemas.ENTRY_WRITER_SCHEMAS.items():
        assert isinstance(schema, pa.Schema), name
