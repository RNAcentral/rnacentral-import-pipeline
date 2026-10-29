from rnacentral_pipeline.databases.helpers import gtdb
from rnacentral_pipeline.databases.mgnify import prepare


def test_marine_catalogue_falls_back_to_marine_metagenome(monkeypatch):
    monkeypatch.setattr(prepare, "get_so_type", lambda _: {"RF00001": "SO:0000652"})
    for getter in ["species", "genus", "family", "order"]:
        monkeypatch.setattr(gtdb, f"get_inferred_{getter}_taxid", lambda _: None)
    data = {
        "data": [
            {
                "inferredPhylogeny": "GTDB:d__Bacteria;p__;c__;o__;f__;g__;s__",
                "additionalAnnotations": {"catalog_name": "marine genome catalogue"},
                "sourceModel": "RFAM:RF00001",
            }
        ],
        "metaData": {},
    }
    entry = prepare.prepare_mgnify_data(data, "")["data"][0]
    assert entry["taxonId"] == "NCBITaxon:408172"
