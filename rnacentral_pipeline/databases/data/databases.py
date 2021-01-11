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
        if name == 'tmrna-website':
            return cls.tmrna_website
        if name == 'snopydb':
            return cls.snopy
        if attribute == 'ensembl/gencode' or attribute == 'ensembl_gencode':
            return cls.gencode
        raise ValueError("Unknown database name %s" % name)

    def normalized(self) -> str:
        if self is Database.gencode:
            return 'ENSEMBL_GENCODE'
        return self.name.upper()
