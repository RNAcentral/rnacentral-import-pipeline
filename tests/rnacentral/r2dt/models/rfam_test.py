# -*- coding: utf-8 -*-

"""
cm-headers.cm is the header of two real Rfam covariance models (RF00002 and
RF00001, from rfam.org/family/<acc>/cm), cut at the "CM" line where
parse_model stops reading. The known info stands in for the rfam_models rows
load_info would read; RF00002's is RNAcentral production's row.
"""

from rnacentral_pipeline.rnacentral.r2dt.data import Source
from rnacentral_pipeline.rnacentral.r2dt.models import rfam

CMS = "data/r2dt/rfam/cm-headers.cm"

# production rfam_models: RF00002 -> so_rna_type SO:0000375, rna_type rRNA
KNOWN = {"RF00002": ("rfam", 131567, "SO:0000375", "rRNA")}


def parse(known, monkeypatch):
    monkeypatch.setattr(rfam, "load_info", lambda db_url: known)
    with open(CMS) as raw:
        return list(rfam.parse(raw, "postgres://unused"))


def test_builds_model_info_from_the_cm_header(monkeypatch):
    (model,) = parse(KNOWN, monkeypatch)
    assert model.model_name == "RF00002"  # ACC, not NAME (5_8S_rRNA)
    assert model.so_rna_type == "SO:0000375"
    assert model.taxid == 131567
    assert model.source == Source.rfam
    # CLEN of the current Rfam CM. Production's r2dt_models row says 155: it
    # was loaded from an earlier Rfam release.
    assert model.length == 154


def test_families_missing_from_rfam_models_are_skipped(monkeypatch):
    # RF00001 is in the file but not in the known info
    assert [m.model_name for m in parse(KNOWN, monkeypatch)] == ["RF00002"]
    assert parse({}, monkeypatch) == []
