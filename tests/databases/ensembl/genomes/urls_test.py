# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import typing as ty
from functools import lru_cache
from io import StringIO

import pytest

from rnacentral_pipeline.databases.ensembl.data import Division, FtpInfo
from rnacentral_pipeline.databases.ensembl.genomes import urls


@lru_cache()
def species(division: Division) -> ty.List[str]:
    found = urls.urls_for(division, "ftp.ensemblgenomes.org")
    return [f[1] for f in found]


@pytest.mark.parametrize(
    "division,expected,found",
    [
        (Division.fungi, "aspergillus_oryzae_gca_002007945", False),
        (Division.fungi, "aspergillus_fumigatus_var_rp_2014_gca_000731615", False),
        (Division.fungi, "aspergillus_oryzae", True),
    ],
)
def test_can_generate_reasonable_species(division, expected, found):
    if found:
        assert expected in species(division)
    else:
        assert expected not in species(division)


def test_can_generate_current_urls_from_species_metadata():
    metadata = StringIO(
        """
[
  {
    "organism": {"name": "aspergillus_oryzae", "url_name": "Aspergillus_oryzae"},
    "assembly": {"assembly_default": "ASM18445v3"},
    "databases": [{"dbname": "aspergillus_oryzae_core_62_1"}]
  }
]
""".strip()
    )

    found = list(
        urls.generate_paths(
            ftp=None,
            division=Division.fungi,
            base="ftp://ftp.ensemblgenomes.org/pub",
            release="current",
            handle=metadata,
        )
    )

    assert found == [
        FtpInfo(
            division=Division.fungi,
            species="aspergillus_oryzae",
            data_files=(
                "ftp://ftp.ensemblgenomes.org/pub/current/fungi/embl/"
                "aspergillus_oryzae/Aspergillus_oryzae.ASM18445v3.62.*.dat.gz"
            ),
            gff_file=(
                "ftp://ftp.ensemblgenomes.org/pub/current/fungi/gff3/"
                "aspergillus_oryzae/Aspergillus_oryzae.ASM18445v3.62.gff3.gz"
            ),
        )
    ]


def test_can_extract_release_suffix_from_database_name():
    entry = {
        "organism": {"name": "aspergillus_oryzae"},
        "databases": [{"dbname": "aspergillus_oryzae_core_62_1"}],
    }

    assert urls.release_suffix(entry) == "62"
