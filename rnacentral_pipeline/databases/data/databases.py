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


@enum.unique
class Database(enum.Enum):
    """
    This is an enum that is used to represent the databases that RNAcentral
    knows about.
    """

    crw = DatabaseValue(0, "CRW")
    dictybase = DatabaseValue(1, "DictyBase")
    ena = DatabaseValue(2, "ENA")
    ensembl = DatabaseValue(3, "Ensembl")
    ensembl_fungi = DatabaseValue(4, "Ensembl Fungi")
    ensembl_metazoa = DatabaseValue(5, "Ensembl Metazoa")
    ensembl_plants = DatabaseValue(6, "Ensembl Plants")
    ensembl_protists = DatabaseValue(7, "Ensembl Protists")
    evlncrnas = DatabaseValue(53, "EVlncRNAs")
    expression_atlas = DatabaseValue(51, "Expression Atlas")
    five_srrnadb = DatabaseValue(8, "5SrRNAdb")
    flybase = DatabaseValue(9, "FlyBase")
    gencode = DatabaseValue(10, "Ensembl/GENCODE")
    genecards = DatabaseValue(11, "GeneCards")
    greengenes = DatabaseValue(12, "Greengenes")
    gtrnadb = DatabaseValue(13, "GtRNAdb")
    hgnc = DatabaseValue(14, "HGNC")
    intact = DatabaseValue(15, "IntAct")
    lncbase = DatabaseValue(16, "LncBase")
    lncbook = DatabaseValue(17, "LncBook")
    lncipedia = DatabaseValue(18, "LNCipedia")
    lncrnadb = DatabaseValue(19, "lncRNAdb")
    malacards = DatabaseValue(20, "MalaCards")
    mgi = DatabaseValue(21, "MGI")
    mirbase = DatabaseValue(22, "miRBase")
    mirgenedb = DatabaseValue(23, "MirGeneDB")
    modomics = DatabaseValue(24, "Modomics")
    noncode = DatabaseValue(25, "NONCODE")
    pdbe = DatabaseValue(26, "PDBe")
    pirbase = DatabaseValue(27, "PirBase")
    plncdb = DatabaseValue(50, "PLncDB")
    pombase = DatabaseValue(28, "PomBase")
    psicquic = DatabaseValue(48, "PSICQUIC")
    rdp = DatabaseValue(29, "RDP")
    refseq = DatabaseValue(30, "RefSeq")
    ribocentre = DatabaseValue(52, "RiboCentre")
    ribovision = DatabaseValue(49, "RiboVision")
    rfam = DatabaseValue(31, "Rfam")
    rgd = DatabaseValue(32, "RGD")
    sgd = DatabaseValue(33, "SGD")
    silva = DatabaseValue(34, "SILVA")
    snodb = DatabaseValue(35, "snoDB")
    snopy = DatabaseValue(36, "snOPY")
    snorna_database = DatabaseValue(37, "snoRNA Database")
    srpdb = DatabaseValue(38, "SRPDB")
    tair = DatabaseValue(39, "TAIR")
    tarbase = DatabaseValue(40, "TarBase")
    tmrna_website = DatabaseValue(41, "tmRNA Website")
    vega = DatabaseValue(42, "VEGA")
    wormbase = DatabaseValue(43, "WormBase")
    zfin = DatabaseValue(44, "Zfin")
    zwd = DatabaseValue(45, "ZWD")

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

    def normalized(self) -> str:
        if self is Database.gencode:
            return "ENSEMBL_GENCODE"
        return self.name.upper().replace(" ", "_")

    def index(self) -> int:
        return self.value.id

    def pretty(self) -> str:
        return self.value.pretty
