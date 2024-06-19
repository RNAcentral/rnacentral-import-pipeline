# -*- coding: utf-8 -*-

from __future__ import annotations

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

import enum
import typing as ty


class DatabaseValue(ty.NamedTuple):
    id: int
    pretty: str
    descr: str

    def matches(self, key: int | str) -> bool:
        if isinstance(key, int):
            return key == self.id
        if isinstance(key, str):
            return (
                key.lower() == self.pretty.lower() or key.lower() == self.descr.lower()
            )
        raise ValueError(f"Cannot match {key}")


@enum.unique
class Database(enum.Enum):
    """
    This is an enum that is used to represent the databases that RNAcentral
    knows about.
    """

    crw = DatabaseValue(0, "CRW", "CRW")
    dictybase = DatabaseValue(1, "DictyBase", "DICTYBASE")
    ena = DatabaseValue(2, "ENA", "ENA")
    ensembl = DatabaseValue(3, "Ensembl", "ENSEMBL")
    ensembl_fungi = DatabaseValue(4, "Ensembl Fungi", "ENSEMBL_FUNGI")
    ensembl_metazoa = DatabaseValue(5, "Ensembl Metazoa", "ENSEMBL_METAZOA")
    ensembl_plants = DatabaseValue(6, "Ensembl Plants", "ENSEMBL_PLANTS")
    ensembl_protists = DatabaseValue(7, "Ensembl Protists", "ENSEMBL_PROTISTS")
    evlncrnas = DatabaseValue(53, "EVlncRNAs", "EVLNCRNAS")
    expression_atlas = DatabaseValue(51, "Expression Atlas", "Expression Atlas")
    five_srrnadb = DatabaseValue(8, "5SrRNAdb", "5SRRNADB")
    flybase = DatabaseValue(9, "FlyBase", "FLYBASE")
    gencode = DatabaseValue(10, "Ensembl/GENCODE", "ENSEMBL_GENCODE")
    genecards = DatabaseValue(11, "GeneCards", "GENECARDS")
    greengenes = DatabaseValue(12, "Greengenes", "GREENGENES")
    gtrnadb = DatabaseValue(13, "GtRNAdb", "GTRNADB")
    hgnc = DatabaseValue(14, "HGNC", "HGNC")
    intact = DatabaseValue(15, "IntAct", "INTACT")
    lncbase = DatabaseValue(16, "LncBase", "LNCBASE")
    lncbook = DatabaseValue(17, "LncBook", "LNCBOOK")
    lncipedia = DatabaseValue(18, "LNCipedia", "LNCIPEDIA")
    lncrnadb = DatabaseValue(19, "lncRNAdb", "LNCRNADB")
    malacards = DatabaseValue(20, "MalaCards", "MALACARDS")
    mgi = DatabaseValue(21, "MGI", "MGI")
    mgnify = DatabaseValue(55, "MGNIFY", "MGNIFY")
    mirbase = DatabaseValue(22, "miRBase", "MIRBASE")
    mirgenedb = DatabaseValue(23, "MirGeneDB", "MIRGENEDB")
    modomics = DatabaseValue(24, "Modomics", "MODOMICS")
    noncode = DatabaseValue(25, "NONCODE", "NONCODE")
    pdbe = DatabaseValue(26, "PDBe", "PDBE")
    pirbase = DatabaseValue(27, "PirBase", "PIRBASE")
    plncdb = DatabaseValue(50, "PLncDB", "PLNCDB")
    pombase = DatabaseValue(28, "PomBase", "POMBASE")
    psicquic = DatabaseValue(48, "PSICQUIC", "PSICQUIC")
    rdp = DatabaseValue(29, "RDP", "RDP")
    refseq = DatabaseValue(30, "RefSeq", "REFSEQ")
    ribocentre = DatabaseValue(52, "RiboCentre", "RIBOCENTRE")
    ribovision = DatabaseValue(49, "RiboVision", "RIBOVISION")
    rfam = DatabaseValue(31, "Rfam", "RFAM")
    rgd = DatabaseValue(32, "RGD", "RGD")
    sgd = DatabaseValue(33, "SGD", "SGD")
    silva = DatabaseValue(34, "SILVA", "SILVA")
    snodb = DatabaseValue(35, "snoDB", "SNODB")
    snopy = DatabaseValue(36, "snOPY", "SNOPY")
    snorna_database = DatabaseValue(37, "snoRNA Database", "SNORNADB")
    srpdb = DatabaseValue(38, "SRPDB", "SRPDB")
    tair = DatabaseValue(39, "TAIR", "TAIR")
    tarbase = DatabaseValue(40, "TarBase", "TARBASE")
    tmrna_website = DatabaseValue(41, "tmRNA Website", "TMRNA_WEB")
    vega = DatabaseValue(42, "VEGA", "VEGA")
    wormbase = DatabaseValue(43, "WormBase", "WORMBASE")
    zfin = DatabaseValue(44, "Zfin", "ZFIN")
    zwd = DatabaseValue(45, "ZWD", "ZWD")

    @classmethod
    def build(cls, name: str) -> Database:
        if isinstance(name, cls):
            return name
        attribute = name.lower().replace(" ", "_")
        if hasattr(cls, attribute):
            return getattr(cls, attribute)
        if name == "pdb":
            return cls.pdbe
        if name == "tmrna-website" or attribute == "tmrna_web":
            return cls.tmrna_website
        if name == "snopydb":
            return cls.snopy
        if attribute == "ensembl/gencode" or attribute == "ensembl_gencode":
            return cls.gencode
        if attribute == "5srrnadb":
            return cls.five_srrnadb
        if attribute == "snornadb":
            return cls.snorna_database
        raise ValueError("Unknown database name %s" % name)

    @classmethod
    def lookup(cls, db: int | str) -> Database | None:
        for value in list(cls):
            if value.value.matches(db):
                return value
        return None

    def normalized(self) -> str:
        if self is Database.gencode:
            return "ENSEMBL_GENCODE"
        return self.name.upper().replace(" ", "_")

    def index(self) -> int:
        return self.value.id

    def pretty(self) -> str:
        return self.value.pretty
