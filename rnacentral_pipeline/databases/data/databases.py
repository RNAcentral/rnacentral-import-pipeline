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

import enum


@enum.unique
class Database(enum.Enum):
    """
    This is an enum that is used to represent the databases that RNAcentral
    knows about.
    """

    crw              = enum.auto()
    dictybase        = enum.auto()
    ena              = enum.auto()
    ensembl          = enum.auto()
    ensembl_fungi    = enum.auto()
    ensembl_metazoa  = enum.auto()
    ensembl_plants   = enum.auto()
    ensembl_protists = enum.auto()
    five_srrnadb     = enum.auto()
    flybase          = enum.auto()
    gencode          = enum.auto()
    genecards        = enum.auto()
    greengenes       = enum.auto()
    gtrnadb          = enum.auto()
    hgnc             = enum.auto()
    intact           = enum.auto()
    lncbase          = enum.auto()
    lncbook          = enum.auto()
    lncipedia        = enum.auto()
    lncrnadb         = enum.auto()
    malacards        = enum.auto()
    mgi              = enum.auto()
    mirbase          = enum.auto()
    mirgenedb        = enum.auto()
    modomics         = enum.auto()
    noncode          = enum.auto()
    pdbe             = enum.auto()
    pirbase          = enum.auto()
    pombase          = enum.auto()
    rdp              = enum.auto()
    refseq           = enum.auto()
    rfam             = enum.auto()
    rgd              = enum.auto()
    sgd              = enum.auto()
    silva            = enum.auto()
    snodb            = enum.auto()
    snopy            = enum.auto()
    snorna_database  = enum.auto()
    srpdb            = enum.auto()
    tair             = enum.auto()
    tarbase          = enum.auto()
    tmrna_website    = enum.auto()
    vega             = enum.auto()
    wormbase         = enum.auto()
    zfin             = enum.auto()
    zwd              = enum.auto()

    @classmethod
    def build(cls, name) -> "Database":
        if isinstance(name, cls):
            return name
        attribute = name.lower().replace(' ', '_')
        if hasattr(cls, attribute):
            return getattr(cls, attribute)
        if name == 'pdb':
            return cls.pdbe
        if name == 'tmrna-website' or attribute == "tmrna_web":
            return cls.tmrna_website
        if name == 'snopydb':
            return cls.snopy
        if attribute == 'ensembl/gencode' or attribute == 'ensembl_gencode':
            return cls.gencode
        if attribute == "5srrnadb":
            return cls.five_srrnadb
        if attribute == "snornadb":
            return cls.snorna_database
        raise ValueError("Unknown database name %s" % name)

    def normalized(self) -> str:
        if self is Database.gencode:
            return 'ENSEMBL_GENCODE'
        return self.name.upper()

    def pretty(self) -> str:
        if self is Database.crw:
            return 'CRW'
        if self is Database.dictybase:
            return 'DictyBase'
        if self is Database.ena:
            return 'ENA'
        if self is Database.ensembl:
            return 'Ensembl'
        if self is Database.ensembl_fungi:
            return 'Ensembl Funig'
        if self is Database.ensembl_metazoa:
            return 'Ensembl Metazoa'
        if self is Database.ensembl_plants:
            return 'Ensembl Plants'
        if self is Database.ensembl_protists:
            return 'Ensembl Protists'
        if self is Database.five_srrnadb:
            return '5SrRNAdb'
        if self is Database.flybase:
            return 'FlyBase'
        if self is Database.gencode:
            return 'Ensembl/GENCODE'
        if self is Database.genecards:
            return 'GeneCards'
        if self is Database.greengenes:
            return 'Greengenes'
        if self is Database.gtrnadb:
            return 'GtRNAdb'
        if self is Database.hgnc:
            return 'HGNC'
        if self is Database.intact:
            return 'IntAct'
        if self is Database.lncbase:
            return 'LncBase'
        if self is Database.lncbook:
            return 'LncBook'
        if self is Database.lncipedia:
            return 'LNCipedia'
        if self is Database.lncrnadb:
            return 'lncRNAdb'
        if self is Database.malacards:
            return 'MalaCards'
        if self is Database.mgi:
            return 'MGI'
        if self is Database.mirbase:
            return 'miRBase'
        if self is Database.mirgenedb:
            return 'MirGeneDB'
        if self is Database.modomics:
            return 'Modomics'
        if self is Database.noncode:
            return 'NONCODCE'
        if self is Database.pdbe:
            return 'PDBe'
        if self is Database.pirbase:
            return 'PirBase'
        if self is Database.pombase:
            return 'PomBase'
        if self is Database.rdp:
            return 'RDP'
        if self is Database.refseq:
            return 'RefSeq'
        if self is Database.rfam:
            return 'Rfam'
        if self is Database.rgd:
            return 'RGD'
        if self is Database.sgd:
            return 'SGD'
        if self is Database.silva:
            return 'SILVA'
        if self is Database.snodb:
            return 'snoDB'
        if self is Database.snopy:
            return 'snOPY'
        if self is Database.snorna_database:
            return 'snoRNA Database'
        if self is Database.srpdb:
            return 'SRPDB'
        if self is Database.tair:
            return 'TAIR'
        if self is Database.tarbase:
            return 'TarBase'
        if self is Database.tmrna_website:
            return 'tmRNA Website'
        if self is Database.vega:
            return 'VEGA'
        if self is Database.wormbase:
            return 'WormBase'
        if self is Database.zfin:
            return 'Zfin'
        if self is Database.zwd:
            return 'ZWD'
        raise ValueError("No pretty name for %s" % self)
