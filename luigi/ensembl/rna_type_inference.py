# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

from rfam import utils as rfutil


LNC_ALIASES = set([
    'sense_intronic',
    'sense_overlapping',
    'macro_lncRNA',
    'bidirectional_promoter_lncRNA',
    'lincRNA',
])

MITO_TYPES = set([
    'Mt_rRNA',
    'Mt_tRNA',
])

NC_ALIASES = set([
    '3prime_overlapping_ncRNA',
    'known_ncrna',
    'non_coding',
])


MGI_TYPES = {
    'Terc': 'telomerase_RNA',
    'Rn7s6': 'SRP_RNA',
    'Nkx2-2os': 'lncRNA',
    # Nbr1-203
    # Vis1-201
    # 1600017P15Rik-201
}


class UnknownRnaTypeException(Exception):
    pass


class RnaTypeInference(object):
    def __init__(self):
        self.rfam_mapping = rfutil.name_to_isnsdc_type()
        self.fallbacks = [
            ('RFAM_trans_name', self.rfam_type),
            ('MGI_trans_name', self.mouse_type),
        ]

    def rfam_type(self, name):
        rfam_name = name.split('.')[0]
        if rfam_name not in self.rfam_mapping and '-' in rfam_name:
            rfam_name, _ = rfam_name.split('-', 1)
        return self.rfam_mapping.get(rfam_name, None)

    def mouse_type(self, name):
        return MGI_TYPES.get(name, None)

    def compute_fallback_rna_type(self, current):
        xrefs = current.xref_data
        for key, func in self.fallbacks:
            if key not in xrefs:
                continue
            values = {func(v) for v in xrefs[key]}
            values.discard(None)
            if len(values) == 1:
                return values.pop()
        raise UnknownRnaTypeException("Could not determine rna_type of %s",
                                      current)

    def infer_rna_type(self, current):
        base_type = current.note[0]
        if base_type in LNC_ALIASES:
            return 'lncRNA'
        if base_type in NC_ALIASES:
            return 'ncRNA'
        if base_type in MITO_TYPES:
            return base_type.replace('Mt_', '')
        if base_type == 'scaRNA':
            return 'snoRNA'
        if base_type == 'misc_RNA':
            return self.compute_fallback_rna_type(current)
        return base_type
