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

import io
import json
import random

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from rnacentral_pipeline.databases.regulondb import parser

GENOME = "".join(random.Random(0).choices("ACGT", k=700))

TUS = [
    "1)id\t2)name\t3)operonId\t4)operonName\t5)tuGenes\t6)promoterId\t7)promoterName"
    "\t8)terminatorIds\t9)tuEvidence\t10)addEvidence\t11)confidenceLevel\t12)PMIDs",
    "TUA\tthrLA\tOPA\tthr\tthrL;thrA\tPA\tthrLp\tTA1,TA2\t\t\tS\t111;222;111",
    "TUB\tspy\tOPB\tspy\tspy\tPB\tspyp\tTB\t\t\tC\t",
    "TUC\tcyaR\tOPC\tcyaR\tcyaR\tPC\tcyaRp\tTC\t\t\tS\t",
    "TUD\tserW\tOPD\tserW\tserW\tPD\tserWp\tTD\t\t\tS\t",
    "TUW\tthrLA\tOPA\tthr\tthrL;thrA\tPA\tthrLp\tTA1\t\t\tW\t",
    "TUL\tleader\tOPL\tleader\t\tPL\tleaderp\tTL\t\t\tS\t",
    "TUR\tryjA\tOPR\tryjA\tryjA\tPR\tryjAp\tTR\t\t\tS\t",
    "TUG\tglyW\tOPG\tglyW\tglyW\tPG\tglyWp\tTG\t\t\tS\t",
]
PROMOTERS = [
    "1)id\t2)name\t3)strand\t4)posTSS",
    "PA\tthrLp\tforward\t31",
    "PB\tspyp\treverse\t340",
    "PC\tcyaRp\tforward\t361",
    "PD\tserWp\treverse\t490",
    "PL\tleaderp\tforward\t520",
    "PR\tryjAp\treverse\t580",
    "PG\tglyWp\tforward\t611",
]
TERMINATORS = [
    "1)id\t2)leftEndPos\t3)rightEndPos\t4)strand",
    "TA1\t191\t200\tforward",
    "TA2\t201\t210\tforward",
    "TB\t231\t240\treverse",
    "TC\t391\t400\tforward",
    "TD\t411\t420\treverse",
    "TL\t541\t550\tforward",
    "TR\t561\t565\treverse",
    "TG\t691\t700\tforward",
]


def feature(kind, start, stop, strand, **qualifiers):
    qualifiers = {k: [v] for k, v in qualifiers.items()}
    return SeqFeature(
        FeatureLocation(start - 1, stop, strand), type=kind, qualifiers=qualifiers
    )


def genome():
    record = SeqRecord(Seq(GENOME), id="U00096.3", name="U00096")
    record.annotations = {
        "molecule_type": "DNA",
        "organism": "Escherichia coli str. K-12 substr. MG1655",
    }
    record.features = [
        feature("source", 1, 700, 1, db_xref="taxon:511145", mol_type="genomic DNA"),
        feature("CDS", 51, 110, 1, gene="thrL", locus_tag="b0001"),
        feature("CDS", 121, 180, 1, gene="thrA", locus_tag="b0002"),
        feature("CDS", 261, 320, -1, gene="spy", locus_tag="b1743"),
        feature(
            "ncRNA", 361, 390, 1, gene="cyaR", locus_tag="b4438", ncRNA_class="other"
        ),
        feature("tRNA", 431, 480, -1, gene="serW", locus_tag="b0883"),
        feature("ncRNA", 557, 580, -1, gene="ryjA", ncRNA_class="other"),
        feature("tRNA", 621, 660, 1, gene="glyW"),
        feature("ncRNA", 661, 690, 1, gene="3'ETS-glyW", ncRNA_class="other"),
        feature("CDS", 665, 690, 1, gene="nextA"),
    ]
    handle = io.StringIO()
    SeqIO.write(record, handle, "embl")
    handle.seek(0)
    return handle


def tables():
    content = {
        "tus": "# RegulonDB TUSet\n" + "\n".join(TUS),
        "promoters": "\n".join(PROMOTERS),
        "terminators": "\n".join(TERMINATORS),
    }
    data = {name: {"content": value} for name, value in content.items()}
    return io.StringIO(json.dumps({"data": data}))


def revcomp(start, stop):
    return str(Seq(GENOME[start - 1 : stop]).reverse_complement())


@pytest.fixture(scope="module")
def entries():
    return {e.accession: e for e in parser.parse(tables(), genome())}


def test_makes_one_transcript_per_tu_and_terminator(entries):
    utrs = [":five_prime_UTR", ":three_prime_UTR"]
    mrnas = ["TUA_TA1", "TUA_TA2", "TUB_TB"]
    assert set(entries) == {f"REGULONDB:{m}{s}" for m in mrnas for s in [""] + utrs} | {
        "REGULONDB:TUC_TC",
        "REGULONDB:TUD_TD",
        "REGULONDB:TUR_TR",
        "REGULONDB:TUG_TG",
    }


def test_a_forward_mrna_runs_from_its_tss_to_its_terminator(entries):
    mrna = entries["REGULONDB:TUA_TA1"]
    assert mrna.rna_type == "SO:0000234"
    assert mrna.sequence == GENOME[30:200]
    assert mrna.note_data == {"cds": [51, 180]}
    assert mrna.url == "https://regulondb.ccg.unam.mx/operon/OPA"
    assert [r.external_id for r in mrna.references] == ["111", "222"]
    assert mrna.description == "Escherichia coli str. K-12 substr. MG1655 thrLA mRNA"
    assert mrna.ncbi_tax_id == 511145
    assert entries["REGULONDB:TUA_TA1:five_prime_UTR"].sequence == GENOME[30:50]
    assert entries["REGULONDB:TUA_TA1:three_prime_UTR"].sequence == GENOME[180:200]
    assert entries["REGULONDB:TUA_TA2:three_prime_UTR"].sequence == GENOME[180:210]


def test_a_reverse_utr_starts_at_the_tss(entries):
    mrna = entries["REGULONDB:TUB_TB"]
    assert mrna.sequence == revcomp(231, 340)
    assert (mrna.gene, mrna.locus_tag) == ("spy", "b1743")
    region = mrna.regions[0]
    assert (region.assembly_id, region.strand.display_string()) == ("ASM584v2", "-")
    assert [(e.start, e.stop) for e in region.exons] == [(231, 340)]
    assert entries["REGULONDB:TUB_TB:five_prime_UTR"].sequence == revcomp(321, 340)
    assert entries["REGULONDB:TUB_TB:three_prime_UTR"].sequence == revcomp(231, 260)


def test_srna_and_trna_units_keep_their_own_type(entries):
    srna = entries["REGULONDB:TUC_TC"]
    assert (srna.rna_type, srna.sequence) == ("SO:0000655", GENOME[360:400])
    assert srna.related_sequences == []
    trna = entries["REGULONDB:TUD_TD"]
    assert (trna.rna_type, trna.sequence) == ("SO:0000185", revcomp(411, 490))
    assert trna.description.endswith("serW precursor RNA")


def test_an_srna_annotated_past_its_terminator_is_kept(entries):
    srna = entries["REGULONDB:TUR_TR"]
    assert (srna.rna_type, srna.sequence) == ("SO:0000655", revcomp(561, 580))


def test_genes_outside_the_tus_gene_list_are_ignored(entries):
    assert entries["REGULONDB:TUG_TG"].rna_type == "SO:0000185"
