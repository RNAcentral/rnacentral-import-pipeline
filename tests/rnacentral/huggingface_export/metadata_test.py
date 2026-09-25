# -*- coding: utf-8 -*-

from rnacentral_pipeline.rnacentral.huggingface_export import metadata


class FakeCursor:
    """Records the executed query and returns a fixed count, no live database."""

    def __init__(self, count):
        self._count = count
        self.query = None

    def execute(self, query, params=None):
        self.query = query

    def fetchone(self):
        return (self._count,)


def test_get_active_sequence_count_queries_the_urs_column_not_upi():
    # rna and rnc_rna_precomputed's id column is `urs`; it was renamed from
    # `upi` and this query was missed, so it 500'd with "column does not exist".
    cur = FakeCursor(42)

    result = metadata.get_active_sequence_count(cur)

    assert "upi" not in cur.query
    assert "urs" in cur.query
    assert result == 42
