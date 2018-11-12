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

import re
import logging

from rnacentral_pipeline.databases.rfam import utils as rfutil

LOGGER = logging.getLogger(__name__)


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
}


class UnknownRnaTypeException(Exception):
    pass


class RnaTypeInference(object):
    def __init__(self, family_handle):
        self.rfam_mapping = rfutil.name_to_insdc_type(family_handle)
        family_handle.seek(0)
        self.rfam_mapping.update(rfutil.id_to_insdc_type(family_handle))
        self.fallbacks = [
            ('RFAM_trans_name', self.rfam_type),
            ('MGI_trans_name', self.mouse_type),
        ]

    def rfam_xref(self, entry):
        key = self.fallbacks[0][0]
        if key in entry.xref_data:
            return entry.xref_data[key]
        return None

    def rfam_name(self, name):
        if name in self.rfam_mapping:
            return name

        pattern = r'^(.+?)\.\d+-\d+$'
        match = re.match(pattern, name)
        if match:
            return match.group(1)
        return None

    def rfam_type(self, name):
        name = self.rfam_name(name)
        return self.rfam_mapping.get(name, None)

    def mouse_type(self, name):
        pattern = r'^(.+)-.*$'
        match = re.match(pattern, name)
        if name not in MGI_TYPES and match:
            name = match.group(1)
        return MGI_TYPES.get(name, None)

    def compute_fallback_rna_type(self, xrefs):
        for key, func in self.fallbacks:
            if key not in xrefs:
                continue
            values = {func(v) for v in xrefs[key]}
            values.discard(None)
            if len(values) == 1:
                return values.pop()
            if len(values) > 1:
                LOGGER.info("Conflicting rna types for %s according to %s",
                            xrefs, key)

        LOGGER.info("Could not infer an rna-type for %s", xrefs)
        return 'misc_RNA'

    def compute_type(self, xref_data, base_type):
        """
        Determine the RNA type for the given data. This will compute an RNA
        type based upon the base type and if that fails compute a fallback type
        based upon the xref data if it is 'misc_RNA'. if that still fails it
        will use the given type.
        """

        if base_type in LNC_ALIASES:
            return 'lncRNA'
        if base_type in NC_ALIASES:
            return 'other'
        if base_type in MITO_TYPES:
            return base_type.replace('Mt_', '')
        if base_type == 'scaRNA':
            return 'snoRNA'
        if base_type == 'misc_RNA':
            return self.compute_fallback_rna_type(xref_data)
        if base_type == 'miRNA':
            return 'precursor_RNA'
        return base_type

    def correct_spelling(self, rna_type):
        """
        Correct any possible spelling mistakes. Sometimes Ensembl has vaultRNA
        instead of 'vault_RNA' or 'antisense' instead of 'antisense_RNA'. This
        corrects those problems.
        """

        if rna_type.startswith('vault'):
            return 'vault_RNA'
        if rna_type.startswith('antisense'):
            return 'antisense_RNA'
        return rna_type

    def infer_rna_type(self, xref_data, base_type):
        """
        Infer the RNA type for the given xrefs and current base type.
        """

        rna_type = self.compute_type(xref_data, base_type)
        return self.correct_spelling(rna_type)
