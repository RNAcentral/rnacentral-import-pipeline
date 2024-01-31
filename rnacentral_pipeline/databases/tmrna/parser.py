# -*- coding: utf-8 -*-

"""
Copyright [2009-${2024}] EMBL-European Bioinformatics Institute
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
import sys
import typing as ty
from functools import lru_cache

from Bio import Entrez, SeqIO

from rnacentral_pipeline.databases.data import Entry, SequenceFeature

from . import helpers

csv.field_size_limit(sys.maxsize)  # It contains huge rows
Entrez.email = "rnacentral@gmail.com"

# #ID	Form	Segments	Tag	Gencode	Evidence	Taxonomy	Sequence	InstanceCt	Instances	Note
# Paulinella__chromatophora.1	Permuted	Acceptor:1-70,IVS:71-76,Coding:77-275,TagCDS:149-199,CCAequiv:71-73	ANNIVRFSRQAAPVAA*	11	aragorn-1.2.40:108.8,infernal-1.1.2:293.2/cyano_tmRNChromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae	GTTCGGTTATTGCCGAACTAGGTGGCTCACACCAATGTTTCGGACAGCGGTTCGATTCCGCTCAGCTCCAttattaGGGGCTGCAATGGTTTCGACGGGGCATCAGGAGGGTTACTGAAGCCTGCTCGGTAAGAGCAAATTAGTAACAgcgaacaacatcgttcgtttctcccgtcaagcggcccctgtggctgccTGACCCTAGATAGGGAGATGAGGTAAAGTCAGCCTTATAACCCAAATGACTCAAGGGGCCTGTAAGGGCCCCATCATTA	1	CP000815.1/744167-744441
# Paulinella__longichromatophora.1	Permuted	Acceptor:1-70,IVS:71-82,Coding:83-280,TagCDS:155-205,CCAequiv:71-73	ANNIVRFSRQAAPVAA*	11	aragorn-1.2.40:108.8,infernal-1.1.2:279.5/cyano_tmRNA	Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae	GTTCGGTTTTAGCCGAACTAGGTGGCTTACACCAATGTTTCGGACAGCGGTTCGATTCCGCTCAGCTCCActggctagtttcGGGGCTGCAATGGTTTCGACGGGGCATGAGGAAGGTTACTGAAGCCTGCTCGGTAAGAGCAAATTTGTAACAgcgaacaacatcgttcgtttctctcgtcaagctgcccctgtggccgccTGACCTTAGCTAGAGAGATGGGGTAAGTCAGCCTTATAACCCAAATGACTCGTGGGACCTGGAAGGGTCCCTAAGTTT	1	MG264610.1/711476-711755
# Paulinella__micropora.1	Permuted	Acceptor:1-70,IVS:71-76,Coding:77-274,TagCDS:149-199,CCAequiv:71-73	ANNIVRFSRQAALVAA*	11	aragorn-1.2.40:108.5,infernal-1.1.2:286.9/cyano_tmRNA	Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae	GTTCGGTTTTAGCCGAACTAGGTGGCTTACACCAGTGTTTCGGACAGCGGTTCGATTCCGCTCAGCTCCAtcacgaGGGGCTGCAATGGTTTCGACGGGGCATGAGGAAGGTTACTGAAGCCTGCTCGGCAAGAGCAAAATCGTATCTgcgaacaacatcgttcgtttctcccgtcaagctgctcttgtagcagccTGACCTTAGTTAAGGAGATGGGGTAAGTCAGCCTTATAACCCAAATGACTCATGGGACCTGGAAGGGTCCCTAAATTT	2	KY124271.1/713251-713524,LC490351.1/683987-684260
# Paulinella__micropora.2	Permuted	Acceptor:1-70,IVS:71-76,Coding:77-274,TagCDS:149-199,CCAequiv:71-73	ANNIVRFSRQAAPVAA*	11	aragorn-1.2.40:108.5,infernal-1.1.2:288.4/cyano_tmRNA	Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae	GTTCGGTTTTAGCCGAACTAGGTGGCTTACACCAGTGTTTCGGACAGCGGTTCGATTCCGCTCAGCTCCAtcacgaGGGGCTGCAATGGTTTCGACGGGGCATGAGGAAGGTTACTGAAGCCTGCTCGGCAAGAGCAAAATCGTATCTgcgaacaacatcgttcgtttctcccgtcaagctgctcctgtagcagccTGACCTTAGTTAAGGAGATGGGGTAAGTCAGCCTTATAACCCAAATGACTCATGGGACCTGGAAGGGTCCCTAAATTT	2	KX897545.1/713048-713321,MG976688.1/713153-713426
# Schizocladia__ischiensis.1	Standard	Body:1-380,TagCDS:84-128,CCAequiv:381-383	VNNIITFNKKLTFA*	11	aragorn-1.2.40:111.1	Plastids,d__Eukaryota	GGGGCTGTTTTGGTTTTGACATTTAAAATGAAATAAATTAATAAGCAGAATACAATAGACATTGTATCCAATTAAGAATAATTgtaaacaacattattacatttaataaaaaactaacttttgcaTAAAATTTAGGAGTTTTTTATGGTTAATTTAATATAGAATTAACTTATATGATAAAACTATTGCTCTAAAATTTAATACTTTTTAGGTAAGTACAATCAACTATAAAATAATTTACTATTTTTTCCATTTGTTATAAAGATTAAATTAATCTCTGATAAAATTCACTAAAATAAAATCTAAAAATTAACTAAATCTGTGAATTAAAATAATTCATTTTATTTAAATGGACGTGGGTTCAATTCCCACCAGCTCCAata	1	NC_053868.1/437-819
# Cryptomonas__curvata.1	Standard	Body:1-335,TagCDS:95-139,CCAequiv:336-338	ANNILSFERKLALV*	11	aragorn-1.2.40:109.0,infernal-1.1.2:109.0/tmRNA	Plastids,d__Eukaryota,c__Cryptophyceae,o__Cryptomonadales,f__Cryptomonadaceae	GGGACTGTTCAGGTATCGACACTCTCTAAAATTTTGTATTATGATTCAAGTCAAGCTTAAATTTTCTTGTAAAACAAAATTTAAAACTATAAACgcaaacaacattctgtcgtttgaacgcaaacttgctttagtaTAAACCTAAAACTAGTTTAAATTATAAAACACATAAGTCGAATAACAGGAAGTTTCTAAATAACTAAAAACTATTGCAATTCCCGACAATCTGAATGAGATCTAAAAATAGTTGCTTTAAATATTGAATAAAGCTAAACTTGTGAATGAATATATAAACGTTGAGCGAGTGGACGTGGGTTCAATTCCCACCAGTTCCAtat	1	NC_035720.1/52138-51801
# Cryptomonas__paramecium.1	Standard	Body:1-298,TagCDS:88-156,CCAequiv:299-301	ASNIVSFQKSPSLASKLFSHRI*	11	aragorn-1.2.40:101.2	Plastids,d__Eukaryota,c__Cryptophyceae,o__Cryptomonadales,f__Cryptomonadaceae	GGGGCTGTAAGGCGTCGACATTCATGGAACTGAGAAGTCAAACAAGTCAAGTTTAAACCATCTTGTAAAAATGGAGTTTAATTAAATgcaagcaacatagtttcatttcaaaagtccccttcactcgcttctaagttattctcgcatagaattTAGAACTCATTAATTACGAAAGAAAAATTCTCAAGAAATACACCAAAAAGAGAGAAGTCGCCTTATATTTTGAATAAGGCTATACTTGTAACTATAGACGCTAAAAGATAATGAATGGACGTGAGTTCAAATCTCACCAGCTCCAaat	1	NC_013703.1/6162-5862
# Cryptomonas__pyrenoidifera.1	Standard	Body:1-337,TagCDS:95-139,CCAequiv:338-340	ANNILSFERKFALV*	11	aragorn-1.2.40:107.3	Plastids,d__Eukaryota,c__Cryptophyceae,o__Cryptomonadales,f__Cryptomonadaceae	GGGACTGTTCAGGTATCGACACTTTCTATAATTTTATACTATGATTCAAGTCAAGCGTAAACTTTCTTGTAAAACTAAGTTTAAAAATACAAATgcaaacaatattctatcgtttgaacgcaaatttgctttagtaTAAACCTAAAAAAATTATGGCTTAAACTATAACAAATACTAATTAACATAGGAGATTTTTAATTCTTAAAAATCAGTATTATCCAAAATGAATTCGGTGAGATAAAAAAATAGTTGCTTTAAATATGGAATAAAGCTAAACTTGTGAATGAATATATATAAATTAAGCGAGTGGACGTGGGTTCAATTCCCACCAGTTCCAtag	1	NC_069042.1/115811-116150
# Chroomonas__placoidea.1	Standard	Body:1-360,TagCDS:95-139,CCAequiv:361-363	ANNIIPFSRKVALV*	11	aragorn-1.2.40:109.5,infernal-1.1.2:108.6/tmRNA	Plastids,d__Eukaryota,c__Cryptophyceae,o__Pyrenomonadales,f__Chroomonadaceae	GGGGCTGTAAAGGTATCGACACTTTTAACAAAATAATAGTATGATTCAAGTCAAGATCGAGTATATCTTGTAAATAAGGCTCAAAACAATAAATgcaaataatataatacctttctctcgtaaggtagcgttagtaTAGAAATCAGTTTTTATACTTTATAAAAATAGTGCTTGACAAGTAATTAGTTATAATTTTTACTGGGGATAGTAACTACATTACTCTAAAACTAAATATTATATATCCCATCTTGAAACCTAATATGCTAAAAAAAAAGCTTTAATTATTGAATAAAGCTAAACTTGTGAACGAGTATATTATAAAGTTGGGAGTGGACGTGGGTTCAAATCCCACCAGCTCCAaaa	1	NC_035721.1/58418-58056

name_mapping = {
    "Acceptor": "tmrna_acceptor",
    "Body": "tmrna_body",
    "CCAequiv": "tmrna_ccaequiv",
    "Coding": "tmrna_coding_region",
    "Exon1": "tmrna_exon",
    "Exon2": "tmrna_exon",
    "GpI": "tmrna_gpi",
    "IVS": "tmrna_ivs",
    "TagCDS": "tmrna_tagcds",
}

NO_CODING_SEQUENCE = ["frameshift", "undetermined"]

CODING_TAGS = ["TagCDS", "Exon1", "Exon2"]

ORGANELLE = ["Chromatophores", "Mitochondria", "Plastids"]


def cleanup_phylogeny(phylogeny: str) -> str:
    return phylogeny


def inferred_species(id: str) -> str:
    raw = id.split(".")[0]
    return raw.replace("__", " ")


def features(raw: ty.Dict[str, str]) -> ty.List[SequenceFeature]:
    features = []
    segements = raw["Segments"].split(",")
    for segement in segements:
        name, range = segement.split(":")
        start, stop = range.split("-")
        # Convert to zero based
        location = [int(start) - 1, int(stop)]
        metadata = {}
        if name in CODING_TAGS:
            if raw["Tag"] not in NO_CODING_SEQUENCE:
                metadata["orf_summary"] = raw["Tag"]
            else:
                coding = raw["Tag"]
                metadata["orf_summary"] = "Coding sequence"
                metadata["has_stop"] = coding[-1] == "*"
                if coding[-1] == "*":
                    coding = coding[0:-1]
                metadata["coding_sequence"] = coding
        feature = SequenceFeature(
            name=name_mapping[name],
            feature_type=name_mapping[name],
            location=location,
            sequence="",
            provider="TMRNA_WEB",
            metadata=metadata,
        )
        features.append(feature)
    return features


def parse(raw: ty.IO) -> ty.Iterable[Entry]:
    reader = csv.DictReader(raw, delimiter="\t")
    rows = list(reader)
    taxid_mapping = helpers.taxid_mapping(rows)
    for row in rows:
        accessions = row["Instances"].split(",")
        assert len(accessions) == int(row["InstanceCt"])
        species = inferred_species(row["#ID"])
        for accession in accessions:
            accession_parts = accession.split("/")
            tax_id = taxid_mapping[accession]
            inferred = cleanup_phylogeny(row["Taxonomy"])
            note = {
                "tmrna_form": row["Form"],
            }
            yield Entry(
                accession=f"tmrna:{accession}",
                primary_id=f"tmrna:{row['#ID']}",
                ncbi_tax_id=tax_id,
                database="TMRNA_WEB",
                sequence=row["Sequence"].upper(),
                regions=[],
                rna_type="SO:0000584",
                url="",
                seq_version="1",
                note_data=note,
                parent_accession=accession_parts[0],
                features=features(row),
                inference=inferred,
                description=f"{species} {row['Form']} tmRNA",
            )
