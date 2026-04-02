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

    ena = DatabaseValue(1, "ENA", "ENA")
    rfam = DatabaseValue(2, "RFAM", "RFAM")
    srpdb = DatabaseValue(3, "SRPDB", "SRPDB")
    mirbase = DatabaseValue(4, "MIRBASE", "MIRBASE")
    vega = DatabaseValue(5, "VEGA", "VEGA")
    tmrna_website = DatabaseValue(6, "tmRNA Website", "TMRNA_WEB")
    lncrnadb = DatabaseValue(7, "lncRNAdb", "LNCRNADB")
    gtrnadb = DatabaseValue(8, "GtRNAdb", "GTRNADB")
    refseq = DatabaseValue(9, "RefSeq", "REFSEQ")
    rdp = DatabaseValue(10, "RDP", "RDP")
    pdbe = DatabaseValue(11, "PDBe", "PDBE")
    snopy = DatabaseValue(12, "snOPY", "SNOPY")
    greengenes = DatabaseValue(13, "Greengenes", "GREENGENES")
    tair = DatabaseValue(14, "TAIR", "TAIR")
    wormbase = DatabaseValue(15, "WormBase", "WORMBASE")
    sgd = DatabaseValue(16, "SGD", "SGD")
    silva = DatabaseValue(17, "SILVA", "SILVA")
    pombase = DatabaseValue(18, "PomBase", "POMBASE")
    dictybase = DatabaseValue(19, "dictyBase", "DICTYBASE")
    lncipedia = DatabaseValue(20, "LNCipedia", "LNCIPEDIA")
    noncode = DatabaseValue(21, "NONCODE", "NONCODE")
    modomics = DatabaseValue(22, "Modomics", "MODOMICS")
    hgnc = DatabaseValue(23, "HUGO Gene Nomenclature Committee", "HGNC")
    flybase = DatabaseValue(24, "FlyBase", "FLYBASE")
    ensembl = DatabaseValue(25, "Ensembl", "ENSEMBL")
    mgi = DatabaseValue(27, "Mouse Genome Database", "MGI")
    rgd = DatabaseValue(28, "Rat Genome Database", "RGD")
    tarbase = DatabaseValue(29, "TarBase", "TARBASE")
    zwd = DatabaseValue(30, "ZWD", "ZWD")
    ensembl_plants = DatabaseValue(31, "Ensembl Plants", "ENSEMBL_PLANTS")
    lncbase = DatabaseValue(32, "LncBase", "LNCBASE")
    lncbook = DatabaseValue(33, "LncBook", "LNCBOOK")
    ensembl_metazoa = DatabaseValue(34, "Ensembl Metazoa", "ENSEMBL_METAZOA")
    ensembl_protists = DatabaseValue(35, "Ensembl Protists", "ENSEMBL_PROTISTS")
    ensembl_fungi = DatabaseValue(36, "Ensembl Fungi", "ENSEMBL_FUNGI")
    snodb = DatabaseValue(37, "snoDB", "SNODB")
    five_srrnadb = DatabaseValue(38, "5SrRNAdb", "5SRRNADB")
    mirgenedb = DatabaseValue(39, "MirGeneDB", "MIRGENEDB")
    malacards = DatabaseValue(40, "MalaCards", "MALACARDS")
    genecards = DatabaseValue(41, "GeneCards", "GENECARDS")
    intact = DatabaseValue(42, "IntAct", "INTACT")
    snorna_database = DatabaseValue(43, "snoRNA Database", "SNORNADB")
    zfin = DatabaseValue(44, "ZFIN", "ZFIN")
    crw = DatabaseValue(45, "CRW", "CRW")
    pirbase = DatabaseValue(46, "PirBase", "PIRBASE")
    gencode = DatabaseValue(47, "ENSEMBL_GENCODE", "ENSEMBL_GENCODE")
    psicquic = DatabaseValue(48, "PSICQUIC", "PSICQUIC")
    ribovision = DatabaseValue(49, "RiboVision", "RIBOVISION")
    plncdb = DatabaseValue(50, "PLncDB", "PLNCDB")
    expression_atlas = DatabaseValue(51, "Expression Atlas", "EXPRESSION_ATLAS")
    ribocentre = DatabaseValue(52, "RiboCentre", "RIBOCENTRE")
    evlncrnas = DatabaseValue(53, "EVlncRNAs", "EVLNCRNAS")
    mgnify = DatabaseValue(55, "MGnify", "MGNIFY")
    mirtrondb = DatabaseValue(56, "mirtronDB", "MIRTRONDB")
    japonicusdb = DatabaseValue(58, "JaponicusDB", "JAPONICUSDB")

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
