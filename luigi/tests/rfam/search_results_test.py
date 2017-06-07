import tempfile
import unittest as ut

import pytest

from rfam.search_results import parse

##target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
##------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
#URS000046C624        -         5S_rRNA              RF00001    cm        1      119      194      312      +    no    1 0.55   0.0  131.6   6.9e-29 !   Avena murphyi partial 5S ribosomal RNA



class ParsingTest(ut.TestCase):
    def test_can_parse_valid_file(self):
        data = parse('data/rfam-hits.tbl')
        assert len(data) == 998

    def test_produces_expected_data(self):
        data = parse('data/rfam-hits.tbl')
        assert data[0] == {
            'mdl_from': 1,
            'mdl_to': 119,
            'seq_from': 194,
            'seq_to': 312',
            'strand': '+',
            'urs': 'URS000046C624',
            'seq_acc': None,
            'rfam_acc': 'RF00001',
            'rfam_id': '5S_rRNA',
            'trunc': 'NONE',
            'pass': 1,
            'gc': 0.55,
            'bias': 0.0,
            'score': 131.6,
            'e_value': 6.9e-29,
            'inc': '!',
            'description': 'Avena murphyi partial 5S ribosomal RNA',
        }

    def test_fails_empty_file(self):
        with tempfile.NamedTemporaryFile() as tfile:
            with pytest.raises():
                parse(tfile.name)

    def test_fails_missing_file(self):
        with pytest.raises():
            parse('/something-that/does-not/exists')

    @pytest.mark.skip()
    def test_fails_with_incomplete_file(self):
        pass
