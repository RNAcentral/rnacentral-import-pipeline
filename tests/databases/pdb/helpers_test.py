# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import datetime as dt

import pytest

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.helpers import publications as pubs
from rnacentral_pipeline.databases.pdb import fetch, helpers
from rnacentral_pipeline.databases.pdb.data import ChainInfo


def load(pdb_id: str, chain_id: str) -> ChainInfo:
    chains = fetch.chains({(pdb_id, chain_id)})
    assert len(chains) == 1
    return chains[0]


@pytest.mark.parametrize(
    "product,expected",
    [
        ("SRP RNA DOMAIN IV", "SRP_RNA"),
        ("PRE-TRNA BULGE-HELIX-BULGE MOTIF", "tRNA"),
        ("tmRNA (63-MER)", "tmRNA"),
        ("U2 snRNA", "snRNA"),
        ("HAMMERHEAD RIBOZYME-RNA STRAND", "hammerhead_ribozyme"),
        ("RNA (HEPATITIS DELTA VIRUS GENOMIC RIBOZYME)", "ribozyme"),
        ("U65 H/ACA snoRNA", "ncRNA"),
        ("RNA (5'-R(*CP*GP*CP*GP*AP*AP*UP*UP*AP*GP*CP*G)-3')", "misc_RNA"),
        ("5S RRNA LOOP D/LOOP E", "rRNA"),
        ("5S RIBOSOMAL RNA", "rRNA"),
        ("5.8S RRNA", "rRNA"),
        ("16S RRNA (5'-R(*GP*GP*CP*GP*UP*CP*AP*CP*AP*CP*CP*UP*UP*C)-3')", "rRNA"),
        ("RNA (16S RNA)", "rRNA"),
        ("Helix 39 of 18S rRNA", "rRNA"),
        ("18S ribosomal RNA", "rRNA"),
        ("23S RRNA", "rRNA"),
        ("HELIX 95 OF 23S RRNA", "rRNA"),
        ("23S RRNA FRAGMENT", "rRNA"),
        ("28S rRNA", "rRNA"),
        ("sarcin/ricin 28S rRNA", "rRNA"),
        ("30S 16S ribosomal RNA", "rRNA"),
        ("30S RNA helix 8", "rRNA"),
        ("40S ribosomal RNA fragment", "rRNA"),
        ("40S ribosomal RNA", "rRNA"),
        ("40S WHEAT GERM RIBOSOME protein 4", "rRNA"),
        ("60S ribosomal RNA fragment", "rRNA"),
        ("60S ribosomal RNA", "rRNA"),
    ],
)
def test_can_compute_correct_rna_types(product: str, expected):
    assert helpers.compound_rna_type(product) == expected


@pytest.mark.network
@pytest.mark.parametrize(
    "pdb,chain,expected",
    [
        ("7mky", "A", True),
        ("7lyj", "A", True),
        ("5u3g", "B", True),
        ("2l1v", "A", True),
        ("6var", "A", True),
        ("4y1i", "A", True),
        ("4y1i", "B", True),
        ("4y1j", "A", True),
        ("4y1j", "B", True),
        ("4y1m", "A", True),
        ("4y1m", "B", True),
        ("7mky", "A", True),
        ("7lyj", "A", True),
        ("7mlw", "F", True),
    ],
)
def test_can_detect_if_is_ncrna(pdb, chain, expected):
    info = load(pdb, chain)
    assert helpers.is_ncrna(info) == expected


def test_description_omits_source_when_organism_name_is_missing():
    # organism_scientific_name is None whenever PDBe's API doesn't supply one
    # (e.g. no organism recorded for some viral RNA) - description() used to
    # interpolate that into the string regardless, producing a dangling
    # "... from  (PDB ...)" (or literal "... from None ...") in the stored
    # description. The whole "from <source>" phrase is dropped instead.
    info = ChainInfo(
        pdb_id="3t4b",
        chain_id="A",
        release_date=dt.datetime(2011, 10, 12),
        experimental_method="X-ray diffraction",
        entity_id=1,
        taxids=[32630],
        resolution=3.55,
        title="Crystal Structure of the HCV IRES pseudoknot domain",
        sequence="CCUCCCGGG",
        molecule_names=["HCV IRES pseudoknot domain plus crystallization module"],
        molecule_type="RNA",
        organism_scientific_name=None,
    )
    assert (
        helpers.description(info)
        == "HCV IRES pseudoknot domain plus crystallization module (PDB 3T4B, chain A)"
    )


def publication(**overrides):
    """A PDBe publication record, as fetch.py hands it to as_reference."""
    row = {
        "pubmed_id": None,
        "doi": None,
        "title": "Structures of the bacterial ribosome in classical and hybrid states",
        "journal_info": {"pdb_abbreviation": "Science"},
        "author_list": [],
    }
    row.update(overrides)
    return row


def test_a_pubmed_id_wins_over_a_doi():
    # Production: 4V9D is referenced as PMID 21596992 / doi 10.1126/science.1202692
    ref = helpers.as_reference(
        publication(pubmed_id="21596992", doi="10.1126/science.1202692")
    )
    assert ref == pubs.reference(21596992)


@pytest.mark.parametrize(
    "doi", ["10.1126/science.1202692", "DOI:10.1126/science.1202692"]
)
def test_a_doi_is_used_when_there_is_no_pubmed_id(doi):
    assert helpers.as_reference(publication(doi=doi)) == pubs.reference(
        "doi:10.1126/science.1202692"
    )


def test_without_ids_the_reference_is_built_from_the_record():
    # Older production rows (e.g. reference 446176) keep the commas inside each
    # name; the current code drops them, as its comment intends.
    ref = helpers.as_reference(
        publication(
            author_list=[
                {"full_name": "Zhang, D.", "last_name": "Zhang", "first_name": "D."},
                {"full_name": None, "last_name": "Yan", "first_name": "K."},
                {"full_name": None, "last_name": "Liu", "first_name": None},
            ],
            journal_info={"pdb_abbreviation": "Nat.Struct.Mol.Biol."},
        )
    )
    assert ref == Reference(
        authors="Zhang D., Yan K., Liu",
        location="Nat.Struct.Mol.Biol.",
        title="Structures of the bacterial ribosome in classical and hybrid states",
        pmid=None,
        doi=None,
    )
