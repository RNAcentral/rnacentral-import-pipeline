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

    dictybase      = 0
    ena            = 1
    ensembl        = 2
    ensembl_plants = 3
    flybase        = 4
    gencode        = 5
    greengenes     = 5
    gtrnadb        = 6
    hgnc           = 7
    lncbase        = 8
    lncbook        = 9
    lncipedia      = 10
    lncrnadb       = 11
    mgi            = 12
    mirbase        = 13
    modomics       = 14
    noncode        = 15
    pdbe           = 16
    pombase        = 17
    rdp            = 18
    refseq         = 19
    rfam           = 20
    rgd            = 21
    sgd            = 22
    silva          = 23
    snopy          = 24
    srpdb          = 25
    tair           = 26
    tarbase        = 27
    tmrna_website  = 28
    vega           = 29
    wormbase       = 30
    zwd            = 31

    @classmethod 
    def build(cls, name):
        if isinstance(value, cls):
            return value
        attribute = name.lower().replace(' ', '_')
        if hasattr(cls, attribute):
            return getattr(cls, attribute)
        if name == 'pdb':
            return cls.pdbe
        if name == 'tmrna-website':
            return cls.tmrna_website
        if name == 'snopydb':
            return cls.snopy
        raise ValueError("Unknown database name %s" % name)

    def normalized(self):
        return self.name.upper()

    def public_name(self):
        if self is Database.pdbe:
            return 'PDB'
        return self.name.upper()
