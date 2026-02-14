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

import attr
import pytest
import asyncio
from unittest.mock import patch

from rnacentral_pipeline.databases.europepmc import fetch

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.data import IdReference

# Apply epmc marker to all tests in this module
pytestmark = pytest.mark.epmc

# Publication metadata extracted from Europe PMC API to avoid network calls
# Maps PMID to full publication summary data from Europe PMC
PUBLICATIONS = {
    "18372920": {
        "id": "18372920",
        "source": "MED",
        "pmid": "18372920",
        "pmcid": "PMC11968769",
        "doi": "10.1038/onc.2008.72",
        "title": "MicroRNA-21 promotes cell transformation by targeting the programmed cell death 4 gene",
        "authorString": "Lu Z, Liu M, Stribinskis V, Klinge CM, Ramos KS, Colburn NH, Li Y.",
        "journalTitle": "Oncogene",
        "issue": "31",
        "journalVolume": "27",
        "journalIssn": "0950-9232; 1476-5594; ",
        "pubYear": "2008",
        "pageInfo": "4373-4379",
        "pubType": "research support, non-u.s. gov't; research-article; journal article; research support, n.i.h., extramural",
        "isOpenAccess": "N",
        "inEPMC": "Y",
        "inPMC": "Y",
        "hasPDF": "Y",
        "hasBook": "N",
        "hasSuppl": "Y",
        "citedByCount": 514,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2008-03-31",
        "firstIndexDate": "2008-11-06",
        "fullTextIdList": {"fullTextId": ["PMC11968769"]},
    },
    "1903816": {
        "id": "1903816",
        "source": "MED",
        "pmid": "1903816",
        "doi": "10.1016/0022-2836(91)90560-s",
        "title": "Informational redundancy of tRNA(4Ser) and tRNA(7Ser) genes in Drosophila melanogaster and evidence for intergenic recombination",
        "authorString": "Leung J, Sinclair DA, Hayashi S, Tener GM, Grigliatti TA.",
        "journalTitle": "J Mol Biol",
        "issue": "2",
        "journalVolume": "219",
        "journalIssn": "0022-2836; 1089-8638; ",
        "pubYear": "1991",
        "pageInfo": "175-188",
        "pubType": "research support, non-u.s. gov't; journal article",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 7,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "Y",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "1991-05-01",
        "firstIndexDate": "2008-12-24",
        "dbCrossReferenceList": {"dbName": ["EMBL"]},
    },
    "19546886": {
        "id": "19546886",
        "source": "MED",
        "pmid": "19546886",
        "doi": "10.1038/cr.2009.72",
        "title": "Regulation of the cell cycle gene, BTG2, by miR-21 in human laryngeal carcinoma",
        "authorString": "Liu M, Wu H, Liu T, Li Y, Wang F, Wan H, Li X, Tang H.",
        "journalTitle": "Cell Res",
        "issue": "7",
        "journalVolume": "19",
        "journalIssn": "1001-0602; 1748-7838; ",
        "pubYear": "2009",
        "pageInfo": "828-837",
        "pubType": "research support, non-u.s. gov't; journal article",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 140,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2009-07-01",
        "firstIndexDate": "2009-10-06",
    },
    "20533548": {
        "id": "20533548",
        "source": "MED",
        "pmid": "20533548",
        "doi": "10.1002/ijc.25506",
        "title": "Micro-RNA-21 regulates TGF-β-induced myofibroblast differentiation by targeting PDCD4 in tumor-stroma interaction",
        "authorString": "Yao Q, Cao S, Li C, Mengesha A, Kong B, Wei M.",
        "journalTitle": "Int J Cancer",
        "issue": "8",
        "journalVolume": "128",
        "journalIssn": "0020-7136; 1097-0215; ",
        "pubYear": "2011",
        "pageInfo": "1783-1792",
        "pubType": "research support, non-u.s. gov't; journal article",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 99,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2011-04-01",
        "firstIndexDate": "2010-06-10",
    },
    "22034194": {
        "id": "22034194",
        "source": "MED",
        "pmid": "22034194",
        "doi": "10.1002/jcp.24006",
        "title": "MicroRNA-21 represses human cystathionine gamma-lyase expression by targeting at specificity protein-1 in smooth muscle cells",
        "authorString": "Yang G, Pei Y, Cao Q, Wang R.",
        "journalTitle": "J Cell Physiol",
        "issue": "9",
        "journalVolume": "227",
        "journalIssn": "0021-9541; 1097-4652; ",
        "pubYear": "2012",
        "pageInfo": "3192-3200",
        "pubType": "research support, non-u.s. gov't; journal article",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 57,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2012-09-01",
        "firstIndexDate": "2011-10-29",
    },
    "23239100": {
        "id": "23239100",
        "source": "MED",
        "pmid": "23239100",
        "doi": "10.1002/jcb.24479",
        "title": "miR-21 modulates the ERK-MAPK signaling pathway by regulating SPRY2 expression during human mesenchymal stem cell differentiation",
        "authorString": "Mei Y, Bian C, Li J, Du Z, Zhou H, Yang Z, Zhao RC.",
        "journalTitle": "J Cell Biochem",
        "issue": "6",
        "journalVolume": "114",
        "journalIssn": "0730-2312; 1097-4644; ",
        "pubYear": "2013",
        "pageInfo": "1374-1384",
        "pubType": "journal article",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 120,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2013-06-01",
        "firstIndexDate": "2012-12-17",
    },
    "26184978": {
        "id": "26184978",
        "source": "MED",
        "pmid": "26184978",
        "pmcid": "PMC4505325",
        "doi": "10.1038/srep12276",
        "title": "MiR-135b-5p and MiR-499a-3p Promote Cell Proliferation and Migration in Atherosclerosis by Directly Targeting MEF2C",
        "authorString": "Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
        "journalTitle": "Sci Rep",
        "journalVolume": "5",
        "journalIssn": "2045-2322",
        "pubYear": "2015",
        "pageInfo": "12276",
        "pubType": "research support, non-u.s. gov't; research-article; journal article",
        "isOpenAccess": "Y",
        "inEPMC": "Y",
        "inPMC": "Y",
        "hasPDF": "Y",
        "hasBook": "N",
        "hasSuppl": "Y",
        "citedByCount": 76,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "Y",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2015-07-17",
        "firstIndexDate": "2015-07-18",
        "fullTextIdList": {"fullTextId": ["PMC4505325"]},
        "dbCrossReferenceList": {"dbName": ["EMBL", "INTACT"]},
    },
    "27334534": {
        "id": "27334534",
        "source": "MED",
        "pmid": "27334534",
        "doi": "10.1099/ijsem.0.001195",
        "title": "Agathobaculum butyriciproducens gen. nov. \xa0sp. nov., a strict anaerobic, butyrate-producing gut bacterium isolated from human faeces and reclassification of Eubacterium desmolans as Agathobaculum desmolans comb. nov",
        "authorString": "Ahn S, Jin TE, Chang DH, Rhee MS, Kim HJ, Lee SJ, Park DS, Kim BC.",
        "journalTitle": "Int J Syst Evol Microbiol",
        "issue": "9",
        "journalVolume": "66",
        "journalIssn": "1466-5026; 1466-5034; ",
        "pubYear": "2016",
        "pageInfo": "3656-3661",
        "pubType": "journal article",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 35,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "Y",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2016-06-21",
        "firstIndexDate": "2016-06-24",
        "dbCrossReferenceList": {"dbName": ["EMBL"]},
    },
    "27389411": {
        "id": "27389411",
        "source": "MED",
        "pmid": "27389411",
        "pmcid": "PMC5079273",
        "doi": "10.1093/cvr/cvw177",
        "title": "MicroRNA-153 targeting of KCNQ4 contributes to vascular dysfunction in hypertension",
        "authorString": "Carr G, Barrese V, Stott JB, Povstyan OV, Jepps TA, Figueiredo HB, Zheng D, Jamshidi Y, Greenwood IA.",
        "journalTitle": "Cardiovasc Res",
        "issue": "2",
        "journalVolume": "112",
        "journalIssn": "0008-6363; 1755-3245; ",
        "pubYear": "2016",
        "pageInfo": "581-589",
        "pubType": "research-article; journal article",
        "isOpenAccess": "Y",
        "inEPMC": "Y",
        "inPMC": "Y",
        "hasPDF": "Y",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 35,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2016-11-01",
        "firstIndexDate": "2016-07-09",
        "fullTextIdList": {"fullTextId": ["PMC5079273"]},
    },
    "27858507": {
        "id": "27858507",
        "source": "MED",
        "pmid": "27858507",
        "pmcid": "PMC5785218",
        "doi": "10.1080/15476286.2016.1251002",
        "title": "Numerous small hammerhead ribozyme variants associated with Penelope-like retrotransposons cleave RNA as dimers",
        "authorString": "Lünse CE, Weinberg Z, Weinberg Z, Breaker RR.",
        "journalTitle": "RNA Biol",
        "issue": "11",
        "journalVolume": "14",
        "journalIssn": "1547-6286; 1555-8584; ",
        "pubYear": "2017",
        "pageInfo": "1499-1507",
        "pubType": "research support, non-u.s. gov't; research-article; journal article; research support, n.i.h., extramural",
        "isOpenAccess": "N",
        "inEPMC": "Y",
        "inPMC": "Y",
        "hasPDF": "Y",
        "hasBook": "N",
        "hasSuppl": "Y",
        "citedByCount": 15,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "Y",
        "firstPublicationDate": "2017-11-03",
        "firstIndexDate": "2016-11-19",
        "fullTextIdList": {"fullTextId": ["PMC5785218"]},
        "tmAccessionTypeList": {"accessionType": ["pfam"]},
    },
    "28815543": {
        "id": "28815543",
        "source": "MED",
        "pmid": "28815543",
        "pmcid": "PMC5890441",
        "doi": "10.1007/978-981-10-5203-3_9",
        "title": "Understanding the Role of lncRNAs in Nervous System Development.",
        "authorString": "Clark BS, Blackshaw S.",
        "journalTitle": "Adv Exp Med Biol",
        "journalVolume": "1008",
        "journalIssn": "0065-2598; 2214-8019; ",
        "pubYear": "2017",
        "pageInfo": "253-282",
        "pubType": "research-article; review; journal article; research support, n.i.h., extramural",
        "isOpenAccess": "N",
        "inEPMC": "Y",
        "inPMC": "Y",
        "hasPDF": "Y",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 40,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2017-01-01",
        "firstIndexDate": "2017-08-18",
        "fullTextIdList": {"fullTextId": ["PMC5890441"]},
    },
    "375006": {
        "id": "375006",
        "source": "MED",
        "pmid": "375006",
        "doi": "10.1007/bf00271669",
        "title": "Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae",
        "authorString": "Macino G, Tzagoloff A.",
        "journalTitle": "Mol Gen Genet",
        "issue": "2",
        "journalVolume": "169",
        "journalIssn": "0026-8925",
        "pubYear": "1979",
        "pageInfo": "183-188",
        "pubType": "research support, u.s. gov't, non-p.h.s.; journal article",
        "isOpenAccess": "N",
        "inEPMC": "N",
        "inPMC": "N",
        "hasPDF": "N",
        "hasBook": "N",
        "hasSuppl": "N",
        "citedByCount": 20,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "1979-01-01",
        "firstIndexDate": "2008-01-29",
    },
}


@pytest.fixture(scope="module")
def mock_europepmc():
    """Mock Europe PMC API to avoid network calls"""

    async def mock_get_data(id_reference):
        """Mock the async get_data function that fetches from Europe PMC API"""
        # The Europe PMC API returns different results based on the ID type
        # We need to find the publication by PMID, regardless of query type

        # Build the query string to find the PMID
        # IdReference can be PMID, DOI, or PMCID - we need to find matching publication
        external_id = id_reference.external_id
        namespace = id_reference.namespace.name

        # Try to find the publication by matching against different ID types
        matched_pub = None
        for pmid, pub_data in PUBLICATIONS.items():
            if namespace == "pmid" and pub_data.get("pmid") == external_id:
                matched_pub = pub_data
                break
            elif namespace == "doi" and pub_data.get("doi") == external_id:
                matched_pub = pub_data
                break
            elif namespace == "pmcid":
                # PMCID can be with or without "PMC" prefix
                pub_pmcid = pub_data.get("pmcid", "")
                if pub_pmcid == external_id or pub_pmcid == f"PMC{external_id}":
                    matched_pub = pub_data
                    break

        # Return Europe PMC API response structure
        if not matched_pub:
            # No match found - return empty result
            return {
                "hitCount": 0,
                "resultList": {"result": []},
            }

        # Return successful match in Europe PMC API format
        return {
            "hitCount": 1,
            "resultList": {"result": [matched_pub]},
        }

    with patch("rnacentral_pipeline.databases.europepmc.fetch.get_data", side_effect=mock_get_data):
        yield


def lookup(ref_id):
    return attr.asdict(fetch.lookup(IdReference.build(ref_id)))

@pytest.mark.parametrize(
    "raw_id", [28815543, "PMC5890441", "doi:10.1007/978-981-10-5203-3_9", "28815543"]
)
def test_can_fetch_publication(mock_europepmc, raw_id):
    idr = IdReference.build(raw_id)
    res = fetch.summary(idr)
    
    assert res == {
        "id": "28815543",
        "source": "MED",
        "pmid": "28815543",
        "doi": "10.1007/978-981-10-5203-3_9",
        "title": "Understanding the Role of lncRNAs in Nervous System Development.",
        "authorString": "Clark BS, Blackshaw S.",
        "journalTitle": "Adv Exp Med Biol",
        "journalVolume": "1008",
        "pubYear": "2017",
        "journalIssn": "0065-2598; 2214-8019; ",
        "pageInfo": "253-282",
        "pubType": "research-article; review; journal article; research support, n.i.h., extramural",
        "isOpenAccess": "N",
        "inEPMC": "Y",
        "inPMC": "Y",
        "hasPDF": "Y",
        "hasBook": "N",
        "citedByCount": 40,
        "hasReferences": "Y",
        "hasTextMinedTerms": "Y",
        "hasDbCrossReferences": "N",
        "hasLabsLinks": "Y",
        "hasTMAccessionNumbers": "N",
        "firstPublicationDate": "2017-01-01",
        "hasSuppl": "N",
        "pmcid": "PMC5890441",
        "fullTextIdList": {"fullTextId": ["PMC5890441"]},
        "firstIndexDate": "2017-08-18",
    }
def test_complains_given_bad_pmid(mock_europepmc):
    with pytest.raises(Exception):
        fetch.summary(IdReference.build(-1))

@pytest.mark.parametrize(
    "raw_id",
    [27858507, "doi:10.1080/15476286.2016.1251002", "PMCID:PMC5785218"],
)
def test_can_build_reference(mock_europepmc, raw_id):
    assert lookup(raw_id) == attr.asdict(
        Reference(
            authors=u"Lünse CE, Weinberg Z, Weinberg Z, Breaker RR.",
            location="RNA Biol 14(11):1499-1507 (2017)",
            title=(
                "Numerous small hammerhead ribozyme variants associated with "
                "Penelope-like retrotransposons cleave RNA as dimers"
            ),
            pmid=27858507,
            doi="10.1080/15476286.2016.1251002",
            pmcid="PMC5785218",
        )
    )
@pytest.mark.parametrize(
    "pmid,title",
    [
        (
            18372920,
            "MicroRNA-21 promotes cell transformation by targeting the programmed cell death 4 gene",
        ),
        (
            19546886,
            "Regulation of the cell cycle gene, BTG2, by miR-21 in human laryngeal carcinoma",
        ),
        (
            20533548,
            "Micro-RNA-21 regulates TGF-β-induced myofibroblast differentiation by targeting PDCD4 in tumor-stroma interaction",
        ),
        (
            22034194,
            "MicroRNA-21 represses human cystathionine gamma-lyase expression by targeting at specificity protein-1 in smooth muscle cells",
        ),
        (
            23239100,
            "miR-21 modulates the ERK-MAPK signaling pathway by regulating SPRY2 expression during human mesenchymal stem cell differentiation",
        ),
    ],
)
def test_can_deal_with_weird_issues(mock_europepmc, pmid, title):
    data = fetch.lookup(IdReference.build(pmid))
    assert data.title == title

def test_can_deal_with_unicode(mock_europepmc):
    data = fetch.lookup(IdReference.build(27334534))
    assert attr.asdict(data) == {
        "authors": "Ahn S, Jin TE, Chang DH, Rhee MS, Kim HJ, Lee SJ, Park DS, Kim BC.",
        "location": "Int J Syst Evol Microbiol 66(9):3656-3661 (2016)",
        "title": "Agathobaculum butyriciproducens gen. nov. \xa0sp. nov., a strict anaerobic, butyrate-producing gut bacterium isolated from human faeces and reclassification of Eubacterium desmolans as Agathobaculum desmolans comb. nov",
        "pmid": 27334534,
        "doi": "10.1099/ijsem.0.001195",
        "pmcid": None,
    }
def test_builds_correction_location(mock_europepmc):
    assert lookup(26184978) == attr.asdict(
        Reference(
            authors="Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
            location="Sci Rep 5:12276 (2015)",
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
def test_can_handle_missing_volume(mock_europepmc):
    assert lookup(27389411) == attr.asdict(
        Reference(
            authors="Carr G, Barrese V, Stott JB, Povstyan OV, Jepps TA, Figueiredo HB, Zheng D, Jamshidi Y, Greenwood IA.",
            location="Cardiovasc Res 112(2):581-589 (2016)",
            title="MicroRNA-153 targeting of KCNQ4 contributes to vascular dysfunction in hypertension",
            pmid=27389411,
            doi="10.1093/cvr/cvw177",
            pmcid="PMC5079273",
        )
    )
def test_it_can_find_if_duplicate_ext_ids(mock_europepmc):
    assert lookup(375006) == attr.asdict(
        Reference(
            authors="Macino G, Tzagoloff A.",
            location="Mol Gen Genet 169(2):183-188 (1979)",
            title="Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae",
            pmid=375006,
            doi="10.1007/bf00271669",
            pmcid=None,
        )
    )
def test_can_lookup_by_doi(mock_europepmc):
    assert lookup("doi:10.1007/bf00271669") == attr.asdict(
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
    "ref_id",
    [
        26184978,
        "pmid:26184978",
        "PMID:26184978",
        "doi:10.1038/srep12276",
        "DOI:10.1038/srep12276",
    ],
)
def test_can_handle_several_reference_formats(mock_europepmc, ref_id):
    assert lookup(ref_id) == attr.asdict(
        Reference(
            authors="Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
            location="Sci Rep 5:12276 (2015)",
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

def test_caching_works_as_expected(mock_europepmc):
    fetch.summary.cache_clear()
    assert fetch.summary.cache_info().hits == 0
    assert fetch.summary.cache_info().misses == 0
    assert lookup("PMID:375006")
    assert fetch.summary.cache_info().hits == 0
    assert fetch.summary.cache_info().misses == 1
    for count in range(10):
        assert lookup("PMID:375006")
        assert fetch.summary.cache_info().hits == count + 1
        assert fetch.summary.cache_info().misses == 1

def test_can_find_unidexable_publication(mock_europepmc):
    ref = lookup("1903816")
    assert ref == attr.asdict(
        Reference(
            authors="Leung J, Sinclair DA, Hayashi S, Tener GM, Grigliatti TA.",
            location="J Mol Biol 219(2):175-188 (1991)",
            title=(
                "Informational redundancy of tRNA(4Ser) and tRNA(7Ser) genes in "
                "Drosophila melanogaster and evidence for intergenic recombination"
            ),
            pmid=1903816,
            doi="10.1016/0022-2836(91)90560-s",
            pmcid=None,
        )
    )
