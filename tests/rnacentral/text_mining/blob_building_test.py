# -*- coding: utf-8 -*-

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

from xml.etree import cElementTree as ET

import pytest
import textblob as tb

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path


from rnacentral_pipeline.rnacentral.text_mining import blob_building as bb


@pytest.mark.parametrize(
    "xml,result",
    [
        (
            "<sec><p>Alternative Uses Test: pretest and posttest</p></sec>",
            "Alternative Uses Test: pretest and posttest",
        ),
        (
            "<sec><p><bold>Alternative</bold> Uses Test: pretest and posttest</p></sec>",
            "Alternative Uses Test: pretest and posttest",
        ),
        (
            "<sec><p><bold><italic>Alternative Uses Test: pretest and posttest</italic></bold></p>.</sec>",
            "Alternative Uses Test: pretest and posttest",
        ),
        (
            '<sec><p><bold><italic>Alternative Uses Test: pretest and posttest</italic></bold>. A computerized 4-min version of the Alternative Uses Test (AUT; Guilford, <xref rid="B25" ref-type="bibr">1950</xref>, <xref rid="B26" ref-type="bibr">1967</xref>) was administered to measure creative ideation.</p></sec>',
            "Alternative Uses Test: pretest and posttest. A computerized 4-min version of the Alternative Uses Test (AUT; Guilford, 1950, 1967) was administered to measure creative ideation.",
        ),
    ],
)
def test_can_normalize_text(xml, result):
    root = ET.fromstring(xml)
    node = root.find("p")
    assert bb.extract_node_text(node) == result


@pytest.mark.parametrize(
    "xml,result",
    [
        (
            "<sec><p>Alternative Uses Test: pretest and posttest</p></sec>",
            ["Alternative Uses Test: pretest and posttest"],
        ),
        (
            "<sec><p><bold>Alternative</bold> Uses Test: pretest and posttest</p></sec>",
            ["Alternative Uses Test: pretest and posttest"],
        ),
        (
            "<sec><p><bold><italic>Alternative Uses Test: pretest and posttest</italic></bold></p>.</sec>",
            ["Alternative Uses Test: pretest and posttest"],
        ),
        (
            '<sec><p><bold><italic>Alternative Uses Test: pretest and posttest</italic></bold>. A computerized 4-min version of the Alternative Uses Test (AUT; Guilford, <xref rid="B25" ref-type="bibr">1950</xref>, <xref rid="B26" ref-type="bibr">1967</xref>) was administered to measure creative ideation.</p></sec>',
            [
                "Alternative Uses Test: pretest and posttest.",
                "A computerized 4-min version of the Alternative Uses Test (AUT; Guilford, 1950, 1967) was administered to measure creative ideation.",
            ],
        ),
    ],
)
def test_creates_correct_sentances_after_normalizing(xml, result):
    root = ET.fromstring(xml)
    node = root.find("p")
    blob = tb.TextBlob(bb.extract_node_text(node))
    assert list(blob.sentences) == result


@pytest.mark.parametrize(
    "path,expected",
    [
        (
            "data/text-mining/plain-text",
            {
                "pmcid:PMC6186449",
                "pmid:29141248",
                "pmid:30541551",
            },
        ),
        ("data/text-mining/plain-text/PMC6186449.txt", {"pmcid:PMC6186449"}),
        ("data/text-mining/plain-text/PMID:30541551.txt", {"pmid:30541551"}),
        ("data/text-mining/plain-text/PMID:29141248.txt", {"pmid:29141248"}),
        (
            "data/text-mining/full-text/simple.xml",
            {
                "doi:10.3791/51382",
                "doi:10.3389/fnagi.2014.00280",
                "doi:10.3389/fnhum.2014.00827",
                "doi:10.1097/SCS.0000000000000889",
            },
        ),
        ("data/publications/example.xml", {"pmid:26184978"}),
    ],
)
def test_can_build_wrapper_with_expected_ids(path, expected):
    wrapper = bb.build(Path(path))
    ids = {blob.pub_id.normalized_id for blob in wrapper.blobs()}
    assert ids == expected


@pytest.mark.parametrize(
    "path,expected",
    [
        (
            "data/publications/example.xml",
            "MiR-135b-5p and MiR-499a-3p Promote Cell Proliferation and Migration in Atherosclerosis by Directly Targeting MEF2C",
        ),
        (
            "data/publications/",
            "MiR-135b-5p and MiR-499a-3p Promote Cell Proliferation and Migration in Atherosclerosis by Directly Targeting MEF2C",
        ),
        (
            "data/text-mining/plain-text/PMC6186449.txt",
            "It is indicated that the majority of Th2 cytokines,",
        ),
        (
            "data/text-mining/plain-text/PMID:30541551.txt",
            "lncRNA HOTAIR are overexpressed in right CRCs samples (p value < 0 and p value =",
        ),
        (
            "data/text-mining/plain-text/PMID:29141248.txt",
            "hepatocellular carcinoma (HCC) susceptibility in a Chinese population. METHODS:We genotyped three SNPs",
        ),
        (
            "data/text-mining/full-text/simple.xml",
            "is a powerful method to analyze fluorescently labeled cells and",
        ),
        (
            "data/text-mining/full-text/simple.xml",
            "While 3DISCO can be employed on various organs, it has been particularly",
        ),
        (
            "data/text-mining/full-text/simple.xml",
            "tissue clearing techniques including 3DISCO, elegantly demonstrate that high-resolution 3D",
        ),
        (
            "data/text-mining/full-text/simple.xml",
            "Resting-state functional connectivity in anterior cingulate cortex in normal aging",
        ),
        (
            "data/text-mining/full-text/simple.xml",
            "The anterior cingulate cortex (ACC) is an important brain region involved in emotional and cognitive processing.",
        ),
        # ('data/text-mining/full-text/simple.xml',
        #  'differences for resting-state functional connectivity of'),
        (
            "data/text-mining/full-text/simple.xml",
            "BA, Brodmann area; dACC, dorsal ACC; rACC, rostral ACC.",
        ),
        (
            "data/text-mining/full-text/simple.xml",
            "Several limitations should be considered in this study.",
        ),
        # ('data/text-mining/full-text/simple.xml',
        #  'Both versions of the alternative uses task measure creative ideation'),
    ],
)
def test_can_build_wrapper_with_expected_text(path, expected):
    wrapper = bb.build(Path(path))
    found = False
    for blob in wrapper.blobs():
        if expected in blob.blob:
            found = True
    assert found, "Could not find '%s'" % expected
