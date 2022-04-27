# -*- coding: utf-8 -*-

# pylint: disable=no-member,missing-docstring,invalid-name

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from collections import OrderedDict

import pytest
pytestmark = pytest.mark.db

from gffutils import Feature

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import gff3

from .helpers import fetch_coord, fetch_all


def fetch_data(rna_id, assembly):
    data = fetch_coord(rna_id, assembly)
    return [f for f in gff3.regions_as_features(data)]


def fetch_gff_data(assembly):
    data = fetch_all(assembly)
    return list(gff3.regions_as_features(data))


def assert_features_equal(val, ans):
    v = "\n".join(str(v) for v in val)
    a = "\n".join(str(a) for a in ans)
    assert v == a


def test_can_produce_features():
    data = fetch_data("URS000082BE64_9606", "GRCh38")
    ans = [
        Feature(
            seqid="3",
            source="RNAcentral",
            featuretype="transcript",
            start=32676136,
            end=32676226,
            strand="+",
            frame=".",
            attributes=OrderedDict(
                [
                    ("Name", ["URS000082BE64_9606"]),
                    ("type", ["snoRNA"]),
                    ("databases", ["snOPY"]),
                    ("ID", ["URS000082BE64_9606.0"]),
                    ("source", ["alignment"]),
                ]
            ),
        ),
        Feature(
            seqid="3",
            source="RNAcentral",
            featuretype="noncoding_exon",
            start=32676136,
            end=32676226,
            strand="+",
            frame=".",
            attributes=OrderedDict(
                [
                    ("Name", ["URS000082BE64_9606"]),
                    ("type", ["snoRNA"]),
                    ("databases", ["snOPY"]),
                    ("ID", ["URS000082BE64_9606.0:ncRNA_exon1"]),
                    ("Parent", ["URS000082BE64_9606.0"]),
                ]
            ),
        ),
    ]
    assert_features_equal(data, ans)


def test_can_produce_features_with_identity():
    data = fetch_data("URS0000563942_9606", "GRCh38")
    ans = [
        Feature(
            seqid="22",
            source="RNAcentral",
            featuretype="transcript",
            start=40363986,
            end=40364062,
            strand="+",
            frame=".",
            attributes=OrderedDict(
                [
                    ("Name", ["URS0000563942_9606"]),
                    ("type", ["snRNA"]),
                    ("databases", ["ENA"]),
                    ("ID", ["URS0000563942_9606.0"]),
                    ("source", ["alignment"]),
                    # 'identity': ['1.00']  FIXME This should have identity 1.0
                ]
            ),
        ),
        Feature(
            seqid="22",
            source="RNAcentral",
            featuretype="noncoding_exon",
            start=40363986,
            end=40364062,
            strand="+",
            frame=".",
            attributes=OrderedDict(
                [
                    ("Name", ["URS0000563942_9606"]),
                    ("type", ["snRNA"]),
                    ("databases", ["ENA"]),
                    ("ID", ["URS0000563942_9606.0:ncRNA_exon1"]),
                    ("Parent", ["URS0000563942_9606.0"]),
                ]
            ),
        ),
    ]
    assert_features_equal(data, ans)


def test_can_build_feature_for_mapped():
    data = fetch_data("URS0000000098_9606", "GRCh38")
    ans = [
        Feature(
            seqid="10",
            source="RNAcentral",
            featuretype="transcript",
            start=17403508,
            end=17403618,
            strand="+",
            frame=".",
            attributes=OrderedDict(
                [
                    ("Name", ["URS0000000098_9606"]),
                    ("type", ["Y_RNA"]),
                    ("databases", ["ENA"]),
                    ("ID", ["URS0000000098_9606.0"]),
                    ("source", ["alignment"]),
                ]
            ),
        ),
        Feature(
            seqid="10",
            source="RNAcentral",
            featuretype="noncoding_exon",
            start=17403508,
            end=17403618,
            strand="+",
            frame=".",
            attributes=OrderedDict(
                [
                    ("Name", ["URS0000000098_9606"]),
                    ("type", ["Y_RNA"]),
                    ("databases", ["ENA"]),
                    ("ID", ["URS0000000098_9606.0:ncRNA_exon1"]),
                    ("Parent", ["URS0000000098_9606.0"]),
                ]
            ),
        ),
    ]
    assert_features_equal(data, ans)


@pytest.mark.xfail()
def test_can_build_correct_features_for_expert_db_location():
    assert False
