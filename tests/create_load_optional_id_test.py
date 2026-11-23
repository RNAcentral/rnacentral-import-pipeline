import re
from pathlib import Path

import pytest

CREATE_LOAD = (
    Path(__file__).resolve().parents[1] / "files" / "schema" / "create_load.sql"
)


def column_type(table, column):
    sql = CREATE_LOAD.read_text()
    body = re.search(rf"TABLE\s+{table}\s*\((.*?)\);", sql, re.S | re.I).group(1)
    return re.search(rf"^\s*{column}\s+([^,\n]+?)\s+NULL", body, re.M | re.I).group(1)


# circAtlas optional ids list every exon and run past 100 characters
@pytest.mark.parametrize(
    "table", ["load_rnc_accessions", "load_rnacentral", "load_rnacentral_all"]
)
def test_optional_id_is_unbounded(table):
    assert column_type(table, "optional_id").lower() == "text"
