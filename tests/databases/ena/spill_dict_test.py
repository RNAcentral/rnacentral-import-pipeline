# -*- coding: utf-8 -*-

"""
SpillDict holds ena.Context's DR (cross-reference) mapping: in memory until its
estimated payload crosses a threshold, then in a SqliteDict on disk.
"""

import attr

from rnacentral_pipeline.databases.ena.spill_dict import SpillDict, _estimate_size


class Slotted:
    __slots__ = ("database", "primary_id", "count")

    def __init__(self, database, primary_id, count=1):
        self.database = database
        self.primary_id = primary_id
        self.count = count


@attr.s(frozen=True)
class Plain:
    database = attr.ib()
    primary_id = attr.ib()


def test_size_counts_the_key_and_string_attributes():
    assert _estimate_size("AB000001", [Slotted("RFAM", "RF00001", 3)]) == 8 + 4 + 7
    assert _estimate_size("AB000001", [Plain("SILVA", "X1")]) == 8 + 5 + 2
    assert _estimate_size("AB000001", "not a list") == 8


def test_small_mappings_stay_in_memory():
    d = SpillDict(threshold_bytes=1_000)
    d["AB000001"] = [Plain("RFAM", "RF00001")]
    assert not d.spilled
    assert "AB000001" in d and "AB000002" not in d
    assert len(d) == 1
    assert d["AB000001"] == [Plain("RFAM", "RF00001")]


def test_crossing_the_threshold_moves_everything_to_disk():
    d = SpillDict(threshold_bytes=30)
    d["AB000001"] = [Plain("RFAM", "RF00001")]  # 19 bytes: still in memory
    assert not d.spilled
    d["AB000002"] = [Plain("SILVA", "X1")]  # 34 bytes: spills
    assert d.spilled
    assert d._spill_path.exists()
    # everything written before and after the spill is still readable
    d["AB000003"] = [Plain("GTRNADB", "tRNA-Ala")]
    d.commit()
    assert len(d) == 3
    assert d["AB000001"] == [Plain("RFAM", "RF00001")]
    assert "AB000003" in d
    d.close()


def test_close_removes_its_own_temp_file():
    d = SpillDict(threshold_bytes=1)
    d["AB000001"] = [Plain("RFAM", "RF00001")]
    path = d._spill_path
    assert path.exists()
    d.close()
    assert not path.exists()
    d.close()  # idempotent


def test_a_caller_supplied_file_is_kept(tmp_path):
    path = tmp_path / "dr.sqlite"
    with SpillDict(threshold_bytes=1, spill_path=path) as d:
        d["AB000001"] = [Plain("RFAM", "RF00001")]
        assert d.spilled
    assert path.exists()


def test_closing_an_unspilled_dict_is_a_no_op():
    with SpillDict(threshold_bytes=1_000) as d:
        d["AB000001"] = []
    assert not d.spilled
