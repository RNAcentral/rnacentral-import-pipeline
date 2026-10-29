import pytest

from rnacentral_pipeline.databases.helpers import gtdb
from rnacentral_pipeline.databases.helpers import phylogeny as phy

LINEAGE = (
    "d__Bacteria;p__Bacillota;c__Clostridia;o__Peptostreptococcales;"
    "f__Peptostreptococcaceae;g__Romboutsia_E;s__Romboutsia_E sp902362365"
)


@pytest.fixture
def lookups(monkeypatch):
    seen = []

    def taxid(name):
        seen.append(name)
        return 1

    monkeypatch.setattr(phy, "taxid", taxid)
    gtdb.get_inferred_family_taxid.cache_clear()
    gtdb.get_inferred_order_taxid.cache_clear()
    return seen


def test_family_lookup_uses_the_family_name(lookups):
    assert gtdb.get_inferred_family_taxid(LINEAGE) == 1
    assert lookups == ["Peptostreptococcaceae bacterium"]


def test_order_lookup_uses_the_order_name(lookups):
    assert gtdb.get_inferred_order_taxid(LINEAGE) == 1
    assert lookups == ["Peptostreptococcales"]
