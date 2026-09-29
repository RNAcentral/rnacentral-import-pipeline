# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import json

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import mrna

# ENST00000615165 on KI270734.1 (minus strand), as the export query returns it,
# with Ensembl's own UTR coordinates from the GFF3 as the expected answer.
EXONS = [
    (138082, 138667),
    (138743, 138831),
    (142194, 142292),
    (143614, 143789),
    (144749, 144895),
    (145004, 145096),
    (146640, 146721),
    (147624, 147703),
    (148116, 148232),
    (148414, 148478),
    (150350, 150499),
    (150987, 151021),
    (156289, 156497),
    (161689, 161750),
]
RECORD = {
    "transcript": "ENST00000615165.1",
    "rna_id": "URS0000000001_9606",
    "gene": "ENSG00000277196.4",
    "gene_name": None,
    "note": json.dumps({"cds": [138480, 156446], "tags": ["Ensembl_canonical"]}),
    "five_prime_utr": "URS0000000002_9606",
    "three_prime_utr": None,
    "chromosome": "KI270734.1",
    "strand": -1,
    "exons": [{"start": s, "stop": e} for s, e in EXONS],
}


def rows(kind):
    return [
        (f.start, f.end, f.frame)
        for f in mrna.features(RECORD)
        if f.featuretype == kind
    ]


def test_builds_the_mrna_with_its_ids_and_tags():
    [feature] = [f for f in mrna.features(RECORD) if f.featuretype == "mRNA"]
    assert (feature.start, feature.end, feature.strand) == (138082, 161750, "-")
    assert feature.attributes["ID"] == ["ENST00000615165.1"]
    assert feature.attributes["Name"] == ["URS0000000001_9606"]
    assert feature.attributes["tag"] == ["Ensembl_canonical"]


def test_utrs_match_ensembl():
    assert [(s, e) for s, e, _ in rows("five_prime_UTR")] == [
        (156447, 156497),
        (161689, 161750),
    ]
    assert [(s, e) for s, e, _ in rows("three_prime_UTR")] == [(138082, 138479)]


def test_cds_covers_the_rest_with_phases_from_the_start_codon():
    # In transcript order, so on the minus strand the start codon comes first
    cds = rows("CDS")
    assert cds[0] == (156289, 156446, "0")
    assert cds[-1][:2] == (138480, 138667)
    assert sum(e - s + 1 for s, e, _ in cds) % 3 == 0
    done = 0
    for start, stop, phase in cds:
        assert phase == str((3 - done % 3) % 3)
        done += stop - start + 1


def test_every_part_points_at_the_mrna():
    parts = [f for f in mrna.features(RECORD) if f.featuretype != "mRNA"]
    assert {tuple(f.attributes["Parent"]) for f in parts} == {("ENST00000615165.1",)}
    assert len([f for f in parts if f.featuretype == "exon"]) == len(EXONS)


def test_cds_phases_follow_a_start_mid_codon():
    note = {"cds": [138480, 156446], "cds_phase": 1}
    record = dict(RECORD, note=json.dumps(note))
    cds = [
        (f.start, f.end, f.frame)
        for f in mrna.features(record)
        if f.featuretype == "CDS"
    ]
    assert cds[0][2] == "1"
    done = 0
    for start, stop, phase in cds:
        assert phase == str((1 - done) % 3)
        done += stop - start + 1


def test_utr_parts_are_named_by_their_own_urs_when_imported():
    parts = list(mrna.features(RECORD))
    five = [f for f in parts if f.featuretype == "five_prime_UTR"]
    three = [f for f in parts if f.featuretype == "three_prime_UTR"]
    assert {tuple(f.attributes["Name"]) for f in five} == {("URS0000000002_9606",)}
    # A UTR too short to import has no URS to link to
    assert all("Name" not in f.attributes for f in three)
