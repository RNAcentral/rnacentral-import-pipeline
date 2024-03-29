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

import csv
import shutil
import tempfile

import six
import attr
import pytest

from rnacentral_pipeline.databases.data import Reference

from rnacentral_pipeline.databases.europepmc import xml
from rnacentral_pipeline.databases.helpers.publications import reference


@pytest.fixture(scope="module")
def indexed_db_path():
    tmp = tempfile.mkdtemp()
    xml.index_directory("data/publications/", tmp)
    yield tmp
    shutil.rmtree(tmp)


@pytest.fixture(scope="module")
def indexed(indexed_db_path):
    xml.index_directory("data/publications/", indexed_db_path)
    with xml.Cache.build(indexed_db_path).open("r") as db:
        yield db


def test_can_parse_xml_data_correctly():
    with open("data/publications/example.xml", "rb") as raw:
        data = list(xml.parse(raw))
        assert len(data) == 1
        assert attr.asdict(data[0]) == attr.asdict(
            Reference(
                authors="Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
                location="Scientific reports 5:12276 (2015)",
                title=(
                    "MiR-135b-5p and MiR-499a-3p Promote Cell "
                    "Proliferation and Migration in Atherosclerosis by Directly "
                    "Targeting MEF2C"
                ),
                pmid=26184978,
                doi="10.1038/srep12276",
                pmcid="PMC4505325",
            )
        )


@pytest.mark.parametrize(
    "raw_id",
    [
        ("26184978"),
        ("DOI:10.1038/srep12276"),
        ("PMC4505325"),
        ("PMCID:PMC4505325"),
        (26184978),
    ],
)
def test_can_query_indexed_data_correctly(indexed, raw_id):
    id_ref = reference(raw_id)
    assert attr.asdict(indexed.get(id_ref)) == attr.asdict(
        Reference(
            authors="Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
            location="Scientific reports 5:12276 (2015)",
            title=(
                "MiR-135b-5p and MiR-499a-3p Promote Cell "
                "Proliferation and Migration in Atherosclerosis by Directly "
                "Targeting MEF2C"
            ),
            pmid=26184978,
            doi="10.1038/srep12276",
            pmcid="PMC4505325",
        )
    )


@pytest.mark.parametrize(
    "raw_id",
    [
        "doi:10.1007/bf00271669",
        "PMID:375006",
    ],
)
def test_can_query_with_fallback(indexed, raw_id):
    id_ref = reference(raw_id)
    ref = indexed.get(id_ref, allow_fallback=True)
    assert attr.asdict(ref) == attr.asdict(
        Reference(
            authors="Macino G, Tzagoloff A.",
            location="Mol Gen Genet 169(2):183-188 (1979)",
            title="Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae",
            pmid=375006,
            doi="10.1007/bf00271669",
            pmcid=None,
        )
    )


@pytest.mark.parametrize(
    "raw_id",
    [
        26184978,
        "pmid:26184978",
        "PMID:26184978",
        "doi:10.1038/srep12276",
        "DOI:10.1038/srep12276",
    ],
)
def test_produces_expected_writeables(indexed, raw_id):
    id_ref = reference(raw_id)
    ref = indexed.get(id_ref, allow_fallback=True)
    assert list(ref.writeable(["something", "other"])) == [
        [
            "03ea261d3cd947fde1cc8328a4c08127",
            "something",
            "other",
            "Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
            "Scientific reports 5:12276 (2015)",
            (
                "MiR-135b-5p and MiR-499a-3p Promote Cell "
                "Proliferation and Migration in Atherosclerosis by Directly "
                "Targeting MEF2C"
            ),
            26184978,
            "10.1038/srep12276",
        ]
    ]


@pytest.mark.parametrize(
    "raw_id",
    [
        26184978,
        "pmid:26184978",
        "PMID:26184978",
        "doi:10.1038/srep12276",
        "DOI:10.1038/srep12276",
    ],
)
def test_can_write_using_specified_columns(indexed_db_path, raw_id):
    out = six.moves.StringIO()
    content = "something,%s,other\n" % raw_id
    raw = six.moves.StringIO(content)
    xml.write_file_lookup(indexed_db_path, raw, out, column=1, allow_fallback=False)
    out.seek(0)
    assert list(csv.reader(out)) == [
        [
            "03ea261d3cd947fde1cc8328a4c08127",
            "something",
            "other",
            "Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
            "Scientific reports 5:12276 (2015)",
            (
                "MiR-135b-5p and MiR-499a-3p Promote Cell "
                "Proliferation and Migration in Atherosclerosis by Directly "
                "Targeting MEF2C"
            ),
            "26184978",
            "10.1038/srep12276",
        ]
    ]


@pytest.mark.parametrize(
    "raw_id",
    [
        "375006",
        "PMID:375006",
        "DOI:10.1007/bf00271669",
    ],
)
def test_can_write_using_specified_columns_and_allow_fallback(indexed_db_path, raw_id):

    with pytest.raises(xml.UncachedReference):
        with xml.Cache.build(indexed_db_path).open("r") as cache:
            id_ref = reference(raw_id)
            ref = cache.get(id_ref, allow_fallback=False)

    out = six.moves.StringIO()
    raw = six.moves.StringIO("something,%s,other\n" % raw_id)
    xml.write_file_lookup(indexed_db_path, raw, out, column=1, allow_fallback=True)
    out.seek(0)
    assert list(csv.reader(out)) == [
        [
            "da0c9805cab7efd21a0eed0e17a8223a",
            "something",
            "other",
            "Macino G, Tzagoloff A.",
            "Mol Gen Genet 169(2):183-188 (1979)",
            "Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae",
            "375006",
            "10.1007/bf00271669",
        ]
    ]
