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


# Names ENA resolves, with the taxids RNAcentral production's rnc_taxonomy has
# for them. ENA's any-name search ignores case, and so does this stand-in.
PRODUCTION_NAMES = {
    "bacterium": 1869227,
    "bacillota": 1239,
    "clostridia": 186801,
    "peptostreptococcales": 3082720,
    "peptostreptococcaceae": 186804,
    "romboutsia": 1501226,
    "uncultured romboutsia sp.": 1505656,
    "peptostreptococcaceae bacterium": 1904861,
    "euglyphida": 178369,
    "paulinellidae": 39715,
}

TMRNA_LINEAGE = "Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae"


@pytest.fixture
def ena(monkeypatch):
    asked = []

    def taxid(name):
        asked.append(name)
        try:
            return PRODUCTION_NAMES[name.lower()]
        except KeyError:
            raise phy.FailedTaxonId(name)

    monkeypatch.setattr(phy, "taxid", taxid)
    for getter in (
        gtdb.get_by_part,
        gtdb.get_inferred_species_taxid,
        gtdb.get_inferred_genus_taxid,
        gtdb.get_inferred_family_taxid,
        gtdb.get_inferred_order_taxid,
        gtdb.phylogeny_to_taxid,
    ):
        getter.cache_clear()
    return asked


def test_species_falls_back_to_an_uncultured_genus_member(ena):
    # GTDB's placeholder species is not in NCBI; the _E suffix is GTDB's own
    assert gtdb.get_inferred_species_taxid(LINEAGE) == 1505656
    assert ena == ["Romboutsia sp902362365", "uncultured Romboutsia sp."]


def test_genus_lookup_strips_the_gtdb_suffix(ena):
    assert gtdb.get_inferred_genus_taxid(LINEAGE) == 1505656
    assert ena == ["uncultured Romboutsia sp."]


@pytest.mark.parametrize("separator", [",", "; "])
def test_by_part_uses_the_most_specific_known_level(ena, separator):
    lineage = TMRNA_LINEAGE.replace(",", separator)
    assert gtdb.get_by_part(lineage) == 39715  # Paulinellidae


def test_by_part_walks_up_until_a_level_resolves(ena):
    # species and genus are GTDB-only names, the family is in NCBI
    assert gtdb.get_by_part(LINEAGE) == 186804
    assert ena[:2] == ["romboutsia_e sp902362365", "romboutsia_e"]


@pytest.mark.parametrize("fallback,expected", [("bacterium", 1869227), (2, 2)])
def test_by_part_ends_at_the_bacteria_fallback(ena, fallback, expected):
    lineage = "d__Bacteria,p__UBA6262,c__UBA6262,o__WVXT01,f__WVXT01"
    assert gtdb.get_by_part(lineage, bacteria_fallback=fallback) == expected


def test_by_part_gives_up_without_a_known_level(ena):
    assert gtdb.get_by_part("p__UBA6262,c__UBA6262,o__WVXT01") is None


def test_phylogeny_to_taxid_tries_by_part_first(ena):
    # A bacterial lineage always resolves by part (at worst to "bacterium"), so
    # the uncultured-species strategies never run for it here
    assert gtdb.phylogeny_to_taxid(LINEAGE) == 186804


def test_phylogeny_to_taxid_fails_loudly(ena):
    with pytest.raises(ValueError):
        gtdb.phylogeny_to_taxid("p__UBA6262,c__UBA6262,o__WVXT01")
