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
import asyncio
import pytest
from unittest.mock import Mock, patch

from rnacentral_pipeline.databases.ensembl.metadata import karyotypes as karyo


def karyotype(domain, species):
    raw = asyncio.run(karyo.fetch(species, domain))
    return karyo.process(raw)


# Mock data for Ensembl API responses
HOMO_SAPIENS_RESPONSE = {
    "default_coord_system_version": "GRCh38",
    "top_level_region": [
        # MT chromosome (no bands)
        {
            "coord_system": "chromosome",
            "length": 16569,
            "name": "MT"
        },
        # Y chromosome (with detailed cytogenetic bands)
        {
            "coord_system": "chromosome",
            "name": "Y",
            "length": 57227415,
            "bands": [
                {"stain": "acen", "id": "p11.1", "start": 10300001, "end": 10400000},
                {"id": "p11.2", "stain": "gneg", "start": 600001, "end": 10300000},
                {"id": "p11.31", "stain": "gpos50", "start": 300001, "end": 600000},
                {"stain": "gneg", "id": "p11.32", "start": 1, "end": 300000},
                {"start": 10400001, "end": 10600000, "id": "q11.1", "stain": "acen"},
                {"stain": "gneg", "start": 10600001, "end": 12400000, "id": "q11.21"},
                {"id": "q11.221", "end": 17100000, "start": 12400001, "stain": "gpos50"},
                {"end": 19600000, "start": 17100001, "id": "q11.222", "stain": "gneg"},
                {"stain": "gpos50", "id": "q11.223", "end": 23800000, "start": 19600001},
                {"stain": "gneg", "id": "q11.23", "end": 26600000, "start": 23800001},
                {"end": 57227415, "start": 26600001, "stain": "gvar", "id": "q12"},
            ]
        },
        # 192 more regions would be here, but we only need MT and Y for tests
    ] + [{"coord_system": "chromosome", "name": str(i), "length": 1000000} for i in range(1, 193)]
}

GLYCINE_MAX_RESPONSE = {
    "default_coord_system_version": "Glycine_max_v2.1",
    "top_level_region": [
        {
            "coord_system": "chromosome",
            "name": "1",
            "length": 56831624
        }
        # Would have 1189 more chromosomes, but we only need one for the test
    ] + [{"coord_system": "chromosome", "name": str(i), "length": 1000000} for i in range(2, 1191)]
}


@pytest.fixture
def mock_ensembl_api():
    """Mock Ensembl REST API responses for karyotype data"""
    def mock_get(url, **kwargs):
        mock_response = Mock()
        mock_response.raise_for_status = Mock()

        # Match URL to return appropriate response
        if "homo_sapiens" in url:
            mock_response.json.return_value = HOMO_SAPIENS_RESPONSE
        elif "glycine_max" in url:
            mock_response.json.return_value = GLYCINE_MAX_RESPONSE
        else:
            raise ValueError(f"Unexpected URL in test: {url}")

        return mock_response

    with patch("rnacentral_pipeline.databases.ensembl.metadata.karyotypes.requests.get", side_effect=mock_get):
        yield


@pytest.mark.ensembl
def test_builds_empty_karyotype_for_missing_data(mock_ensembl_api):
    _, found = karyotype("ensembl", "glycine_max")
    assert len(found) == 1190
    assert found["1"] == {
        "size": 56831624,
        "bands": [
            {
                "start": 1,
                "end": 56831624,
            }
        ],
    }

@pytest.mark.ensembl
def test_builds_with_known_bands(mock_ensembl_api):
    _, found = karyotype("ensembl", "homo_sapiens")
    assert len(found) == 194
    assert found["MT"] == {
        "size": 16569,
        "bands": [
            {
                "start": 1,
                "end": 16569,
            }
        ],
    }
    assert found["Y"] == {
        "size": 57227415,
        "bands": [
            {"type": "acen", "id": "p11.1", "start": 10300001, "end": 10400000},
            {"id": "p11.2", "end": 10300000, "start": 600001, "type": "gneg"},
            {"end": 600000, "start": 300001, "id": "p11.31", "type": "gpos50"},
            {"end": 300000, "start": 1, "id": "p11.32", "type": "gneg"},
            {"start": 10400001, "end": 10600000, "id": "q11.1", "type": "acen"},
            {"type": "gneg", "start": 10600001, "end": 12400000, "id": "q11.21"},
            {"id": "q11.221", "end": 17100000, "start": 12400001, "type": "gpos50"},
            {"end": 19600000, "start": 17100001, "id": "q11.222", "type": "gneg"},
            {"type": "gpos50", "id": "q11.223", "end": 23800000, "start": 19600001},
            {"type": "gneg", "id": "q11.23", "end": 26600000, "start": 23800001},
            {"type": "gvar", "id": "q12", "start": 26600001, "end": 57227415},
        ],
    }
