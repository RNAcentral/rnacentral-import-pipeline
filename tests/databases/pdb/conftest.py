# -*- coding: utf-8 -*-

"""
Shared test fixtures for PDB tests.

This module provides mock fixtures for PDBe API endpoints to allow tests to run
without network dependencies.
"""

import pytest
from unittest.mock import Mock, patch
from urllib.parse import parse_qs, urlparse


# Mock data for PDBe Search API responses (GET endpoint)
# Structure: https://www.ebi.ac.uk/pdbe/search/pdb/select
MOCK_CHAIN_DATA = {
    "1j5e": {
        "response": {
            "numFound": 1,
            "docs": [
                {
                    "pdb_id": "1j5e",
                    "chain_id": ["A"],
                    "entity_id": 1,
                    "tax_id": [274],
                    "resolution": 3.05,
                    "release_date": "2002-04-12T01:00:00Z",
                    "experimental_method": ["X-ray diffraction"],
                    "title": "Structure of the Thermus thermophilus 30S Ribosomal Subunit",
                    "molecule_sequence": "UUUGUUGGAGAGUUUGAUCCUGGCUCAGGGUGAACGCUGGCGGCGUGCCUAAGACAUGCAAGUCGUGCGGGCCGCGGGGUUUUACUCCGUGGUCAGCGGCGGACGGGUGAGUAACGCGUGGGUUGACCUACCCGGAAGAGGGGGACAACCCGGGGAAACUCGGGCUAAUCCCCCCAUGUGGACCCGCCCCUUGGGGUGUGUGCCAAAGGGCUUUGCCCGCUUCCGGAUGGGGCCGCGUCCCAUCAGCUAGUUGGUGGGGUAAUGGGCCCACCAAGGCGACGACGGGGUAGCCGGUCUGAGAGGAUGGCCGGCCACAGGGGCACUGAGACACGGGCCCCACUCCUACGGGAGGCAGCAGUUAGGAAUCUUCCGCAAUGGGGCGCAAGCCUGACGGAGCGACGCCGCUUGGAGGAAGAAGCCCUUCGGGGUGGUAAACUCCUGAACCCGGGACGAAACCCCCCGACGAGGGGGACUGACGGUACCGGGGGUAAUAGCGCCGGCCAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGCGCGAGCGUUACCCGGAUUCACUGGGCGUAAAGGGCGUGUAGGCGGCCUGGGGGCGUCCCAUGUUGAAAGACCACGGCUCAACCGUGGGGGGAGCGUGGGGAUACGCUCAGGCCUAGACGGUGGGAGAGGGUGGUGGGAAUUCCCGGAGUAGCGGUGAAAUGCGCAGAUACCGGGAGGAACGCCGAUGGCGAAGGCAGCCACCUGGUCCACCCCGUGACGCUGAGGCGCGAAAGCGUGGGGAGCAAACCGGAUUAGAUACCCGGGUAGUUCCACGCCCUAAACGAUUGCGCGCUAGGUCUCUCUGGGUCUCCUGGGGGGCCGAAGCUAACGCGUUAAGCGCGCCGCCUGGGGAGUACGGCCGCAAGGCUGAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGGUUUAAUUCGAAGCAACGCGAAGAACCUUACCAGGCCUUGACAUUGCUAGGGAACCCCGGGUGAAAGCCUGGGGUGCCCCGCGAGGGGGAGCCCUAGCACAGGUGCUGCAUGGCCGUCGUCAGCUCGUGCCGUGAGGUGUUGGGGUUAAGUCCCGCAACGAGCGCAACCCCCCGCCGUUAGUUGCCAGCGGUUCGGCCGGGCACUCUAACGGGACUGCCCGCGAAAGCGGGAGGAAGGAGGGGGACGACGUCUGGUCAGCAUGGCCCCUUACGGCCUGGGCGACACACGUGCUACAAUGGCCACUACAAAGCGAUGCCACCCCGGCAACGGGGAGCUAAUCGCAAAAAGGUUGGGCCCAGUUCGGAUUGGGGGUCUGCAACCCGACCCCAUGAAGCCGGAAUCGCUAGUAAUCGCGGAUCAGCCAUUGCCGCGGUUGAAUACGUUCCCCGGGCCUUGUACACACCGCCCGUCACGCCAUGGGGAGCGGGCCUUACCCGAAGUCGCCGGGAGCCUACGGGCCAAGGGGCCGAGGGUUAGGGCCCGUGACUGGGGGCGAAGUCGUAACAAGGGUAGCUGUACCGGAAGGUGCGGCUGGAUCACCUCCUUUCU",
                    "molecule_name": ["16S ribosomal RNA"],
                    "molecule_type": "RNA",
                    "organism_scientific_name": ["Thermus thermophilus"],
                }
            ],
        }
    },
    "1cq5": {
        "response": {
            "numFound": 1,
            "docs": [
                {
                    "pdb_id": "1cq5",
                    "chain_id": ["A"],
                    "entity_id": 1,
                    "tax_id": [562],
                    "resolution": None,
                    "release_date": "1999-08-23T01:00:00Z",
                    "experimental_method": ["Solution NMR"],
                    "title": "NMR STRUCTURE OF SRP RNA DOMAIN IV",
                    "molecule_sequence": "GGCGUUUACCAGGUCAGGUCCGGAAGGAAGCAGCCAAGGCGCC",
                    "molecule_name": ["SRP RNA DOMAIN IV"],
                    "molecule_type": "RNA",
                    "organism_scientific_name": ["Escherichia coli"],
                }
            ],
        }
    },
    "1s72": {
        "response": {
            "numFound": 31,  # Note: 1S72 has many chains, but we only mock chain 9
            "docs": [
                {
                    "pdb_id": "1s72",
                    "chain_id": ["9"],
                    "entity_id": 2,
                    "tax_id": [2238],
                    "resolution": 2.4,
                    "release_date": "2004-06-15T01:00:00Z",
                    "experimental_method": ["X-ray diffraction"],
                    "title": "REFINED CRYSTAL STRUCTURE OF THE HALOARCULA MARISMORTUI LARGE RIBOSOMAL SUBUNIT AT 2.4 ANGSTROM RESOLUTION",
                    "molecule_sequence": "UUAGGCGGCCACAGCGGUGGGGUUGCCUCCCGUACCCAUCCCGAACACGGAAGAUAAGCCCACCAGCGUUCCGGGGAGUACUGGAGUGCGCGAGCCUCUGGGAAACCCGGUUCGCCGCCACC",
                    "molecule_name": ["5S ribosomal RNA"],
                    "molecule_type": "RNA",
                    "organism_scientific_name": ["Haloarcula marismortui"],
                    "rfam_id": ["5S_rRNA"],
                }
            ],
        }
    },
    "3t4b": {
        "response": {
            "numFound": 1,
            "docs": [
                {
                    "pdb_id": "3t4b",
                    "chain_id": ["A"],
                    "entity_id": 1,
                    "tax_id": [32630, 356418],  # Multiple taxids - testing strange taxid handling
                    "resolution": 3.55,
                    "release_date": "2011-10-12T01:00:00Z",
                    "experimental_method": ["X-ray diffraction"],
                    "title": "Crystal Structure of the HCV IRES pseudoknot domain",
                    "molecule_sequence": "CCUCCCGGGAGAGCCGCUAAGGGGGAAACUCUAUGCGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGA",
                    "molecule_name": ["HCV IRES pseudoknot domain plus crystallization module"],
                    "molecule_type": "RNA",
                    "organism_scientific_name": [],  # Note: No organism name for viral RNA
                }
            ],
        }
    },
}


# Mock data for PDBe Publications API responses (POST endpoint)
# Structure: https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/
MOCK_PUBLICATION_DATA = {
    "1j5e": [
        {
            "pubmed_id": "11014182",
            "doi": "10.1038/35030006",
            "title": "Structure of the 30S ribosomal subunit.",
            "journal_info": {
                "pdb_abbreviation": "Nature",
                "pages": "833-838",
                "volume": "407",
                "year": 2000
            },
            "author_list": [
                {"full_name": "Wimberly B.T.", "last_name": "Wimberly", "initials": "B.T."},
                {"full_name": "Brodersen D.E.", "last_name": "Brodersen", "initials": "D.E."},
                {"full_name": "Clemons Jr. W.M.", "last_name": "Clemons Jr.", "initials": "W.M."},
                {"full_name": "Morgan-Warren R.J.", "last_name": "Morgan-Warren", "initials": "R.J."},
                {"full_name": "Carter A.P.", "last_name": "Carter", "initials": "A.P."},
                {"full_name": "Vonrhein C.", "last_name": "Vonrhein", "initials": "C."},
                {"full_name": "Hartsch T.", "last_name": "Hartsch", "initials": "T."},
                {"full_name": "Ramakrishnan V.", "last_name": "Ramakrishnan", "initials": "V."},
            ],
        },
        {
            "pubmed_id": None,
            "doi": "10.1038/35030019",
            "title": "Functional insights from the structure of the 30S ribosomal subunit and its interactions with antibiotics",
            "journal_info": {
                "pdb_abbreviation": "Nature",
                "pages": "340-348",
                "volume": "407",
                "year": 2000
            },
            "author_list": [
                {"full_name": "Carter A.P.", "last_name": "Carter", "initials": "A.P."},
                {"full_name": "Clemons W.M.", "last_name": "Clemons", "initials": "W.M."},
                {"full_name": "Brodersen D.E.", "last_name": "Brodersen", "initials": "D.E."},
                {"full_name": "Morgan-Warren R.J.", "last_name": "Morgan-Warren", "initials": "R.J."},
                {"full_name": "Wimberly B.T.", "last_name": "Wimberly", "initials": "B.T."},
                {"full_name": "Ramakrishnan V.", "last_name": "Ramakrishnan", "initials": "V."},
            ],
        },
        {
            "pubmed_id": None,
            "doi": None,
            "title": "Structure of a Bacterial 30S Ribosomal Subunit at 5.5 A Resolution",
            "journal_info": {
                "pdb_abbreviation": "Nature",
                "pages": "833-840",
                "volume": "400",
                "year": 1999
            },
            "author_list": [
                {"full_name": "Clemons Jr. W.M.", "last_name": None, "initials": None},
                {"full_name": "May J.L.C.", "last_name": None, "initials": None},
                {"full_name": "Wimberly B.T.", "last_name": None, "initials": None},
                {"full_name": "McCutcheon J.P.", "last_name": None, "initials": None},
                {"full_name": "Capel M.S.", "last_name": None, "initials": None},
                {"full_name": "Ramakrishnan V.", "last_name": None, "initials": None},
            ],
        },
    ],
    "1cq5": [
        {
            "pubmed_id": "10580470",
            "doi": "10.1017/s1355838299991458",
            "title": "Structure of the phylogenetically most conserved domain of SRP RNA.",
            "journal_info": {
                "pdb_abbreviation": "RNA",
                "pages": "1453-1463",
                "volume": "5",
                "year": 1999
            },
            "author_list": [
                {"full_name": "Schmitz U.", "last_name": "Schmitz", "initials": "U."},
                {"full_name": "Behrens S.", "last_name": "Behrens", "initials": "S."},
                {"full_name": "Freymann D.M.", "last_name": "Freymann", "initials": "D.M."},
                {"full_name": "Keenan R.J.", "last_name": "Keenan", "initials": "R.J."},
                {"full_name": "Lukavsky P.", "last_name": "Lukavsky", "initials": "P."},
                {"full_name": "Walter P.", "last_name": "Walter", "initials": "P."},
                {"full_name": "James T.L.", "last_name": "James", "initials": "T.L."},
            ],
        }
    ],
    "1s72": [
        {
            "pubmed_id": "15184028",
            "doi": "10.1016/j.jmb.2004.03.076",
            "title": "The roles of ribosomal proteins in the structure assembly, and evolution of the large ribosomal subunit.",
            "journal_info": {
                "pdb_abbreviation": "J.Mol.Biol.",
                "pages": "1093-1108",
                "volume": "340",
                "year": 2004
            },
            "author_list": [
                {"full_name": "Klein D.J.", "last_name": "Klein", "initials": "D.J."},
                {"full_name": "Moore P.B.", "last_name": "Moore", "initials": "P.B."},
                {"full_name": "Steitz T.A.", "last_name": "Steitz", "initials": "T.A."},
            ],
        }
    ],
    "3t4b": [
        {
            "pubmed_id": "22000514",
            "doi": "10.1016/j.str.2011.08.002",
            "title": "Crystal structure of the HCV IRES central domain reveals strategy for start-codon positioning.",
            "journal_info": {
                "pdb_abbreviation": "Structure",
                "pages": "1456-1466",
                "volume": "19",
                "year": 2011
            },
            "author_list": [
                {"full_name": "Khatter H.", "last_name": "Khatter", "initials": "H."},
                {"full_name": "Myasnikov A.G.", "last_name": "Myasnikov", "initials": "A.G."},
                {"full_name": "Natchiar S.K.", "last_name": "Natchiar", "initials": "S.K."},
                {"full_name": "Klaholz B.P.", "last_name": "Klaholz", "initials": "B.P."},
            ],
        }
    ],
}


# Taxonomy data for core PDB test tax IDs
# Maps tax_id -> ENA API response
MOCK_TAXONOMY_DATA = {
    274: {
        "scientificName": "Thermus thermophilus",
        "lineage": "Bacteria; Deinococcota; Deinococci; Thermales; Thermaceae; Thermus; "
    },
    562: {
        "scientificName": "Escherichia coli",
        "lineage": "Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia; "
    },
    2238: {
        "scientificName": "Haloarcula marismortui",
        "lineage": "Archaea; Euryarchaeota; Stenosarchaea group; Halobacteria; Halobacteriales; Halobacteriaceae; Haloarcula; "
    },
    32630: {
        "scientificName": "Hepatitis C virus genotype 1",
        "lineage": "Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Flasuviricetes; Amarillovirales; Flaviviridae; Hepacivirus; "
    },
    356418: {
        "scientificName": "Hepatitis C virus genotype 1b",
        "lineage": "Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Flasuviricetes; Amarillovirales; Flaviviridae; Hepacivirus; "
    },
}


@pytest.fixture(scope="module")
def mock_pdbe_api():
    """
    Mock PDBe Search and Publications APIs to allow tests to run without network.

    This fixture mocks:
    1. GET requests to https://www.ebi.ac.uk/pdbe/search/pdb/select (chain data)
    2. POST requests to https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/ (publication data)
    3. GET requests to ENA taxonomy API for species lookups

    Only core PDB IDs (1J5E, 1CQ5, 1S72, 3T4B) are fully mocked to support
    non-parametrized tests. Parametrized tests remain network-dependent.
    """

    def mock_get(url, **kwargs):
        """Mock handler for GET requests - handles both PDBe Search and ENA taxonomy APIs"""
        mock_response = Mock()
        mock_response.raise_for_status = Mock()

        # Check if this is an ENA taxonomy API request
        if "ebi.ac.uk/ena/taxonomy/rest" in url or "rest.uniprot.org/taxonomy" in url:
            # Extract tax ID from URL
            # URL format: https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/274
            import re
            tax_id_match = re.search(r'/(\d+)(?:\.json)?$', url)
            if tax_id_match:
                tax_id = int(tax_id_match.group(1))
                if tax_id in MOCK_TAXONOMY_DATA:
                    mock_response.json.return_value = MOCK_TAXONOMY_DATA[tax_id]
                    return mock_response

            # Return empty response if tax ID not mocked
            mock_response.json.return_value = {}
            return mock_response

        # Parse the query parameter to determine which PDB ID is requested
        parsed = urlparse(url)
        query_params = parse_qs(parsed.query)

        if "q" in query_params:
            query = query_params["q"][0]

            # Extract PDB ID from query string (e.g., "pdb_id:1j5e" or "pdb_id:1j5e OR ...")
            # Handle simple single PDB queries
            if "pdb_id:" in query:
                # Extract the first PDB ID from the query
                pdb_id = query.split("pdb_id:")[1].split()[0].lower()

                if pdb_id in MOCK_CHAIN_DATA:
                    mock_response.json.return_value = MOCK_CHAIN_DATA[pdb_id]
                    return mock_response

        # If we get here, the query wasn't recognized - return empty result
        # This allows tests to fail gracefully if they request unmocked PDB IDs
        mock_response.json.return_value = {"response": {"numFound": 0, "docs": []}}
        return mock_response

    def mock_post(url, **kwargs):
        """Mock handler for POST requests to PDBe Publications API"""
        mock_response = Mock()
        mock_response.raise_for_status = Mock()

        # Parse POST data which is comma-separated PDB IDs
        if "data" in kwargs:
            pdb_ids = kwargs["data"].split(",")
            result = {}

            for pdb_id in pdb_ids:
                pdb_id = pdb_id.strip().lower()
                if pdb_id in MOCK_PUBLICATION_DATA:
                    result[pdb_id] = MOCK_PUBLICATION_DATA[pdb_id]

            mock_response.json.return_value = result
            return mock_response

        # Empty result if no data provided
        mock_response.json.return_value = {}
        return mock_response

    # Patch both GET and POST at multiple module levels to catch all API calls
    with patch("rnacentral_pipeline.databases.pdb.fetch.requests.get", side_effect=mock_get), \
         patch("rnacentral_pipeline.databases.pdb.fetch.requests.post", side_effect=mock_post), \
         patch("rnacentral_pipeline.databases.helpers.phylogeny.requests.get", side_effect=mock_get):
        yield
