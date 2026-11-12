# -*- coding: utf-8 -*-

from rnacentral_pipeline.rnacentral.genes.random_forest import data


class Cursor:
    def __init__(self, seen):
        self.seen = seen

    def execute(self, sql, args):
        self.seen.append(sql)

    def fetchall(self):
        return []


class Connection:
    def __init__(self, seen):
        self.seen = seen

    def cursor(self, cursor_factory=None):
        return Cursor(self.seen)

    def close(self):
        pass


# Ensembl coding mRNAs and UTRs overlap the ncRNA transcripts and would reshape their genes.
def test_mrnas_and_utrs_are_left_out_of_gene_building(monkeypatch):
    seen = []
    monkeypatch.setattr(data.pg, "connect", lambda conn_str: Connection(seen))
    data.fetch_data(9606, "postgres://")
    assert "not in ('SO:0000234', 'SO:0000204', 'SO:0000205')" in seen[0]
