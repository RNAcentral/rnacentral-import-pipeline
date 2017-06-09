import tempfile
import unittest as ut

import attr
import pytest

from rfam.utils import tbl_iterator, RfamHit


class ParsingTest(ut.TestCase):
    def test_can_parse_valid_file(self):
        data = list(tbl_iterator('data/rfam-hits.tbl'))
        assert len(data) == 998

    def test_produces_expected_data(self):
        data = list(tbl_iterator('data/rfam-hits.tbl'))
        assert attr.asdict(data[0]) == attr.asdict(RfamHit(
            target_name='URS000046C624',
            seq_acc=None,
            rfam_name='5S_rRNA',
            rfam_acc='RF00001',
            mdl='cm',
            mdl_from=1,
            mdl_to=119,
            seq_from=194,
            seq_to=312,
            strand=1,
            trunc=(False, False),
            infernal_pass=1,
            infernal_gc=0.55,
            bias=0.0,
            score=131.6,
            e_value=6.9e-29,
            inc='unique',
            description='Avena murphyi partial 5S ribosomal RNA',
        ))

    def test_produces_nothing_on_empty_file(self):
        with tempfile.NamedTemporaryFile() as tfile:
            assert list(tbl_iterator(tfile.name)) == []

    def test_fails_missing_file(self):
        with pytest.raises(Exception):
            list(tbl_iterator('/something-that/does-not/exists'))

    @pytest.mark.skip()
    def test_fails_with_incomplete_file(self):
        pass
