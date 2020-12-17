# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import typing as ty

from rnacentral_pipeline.rnacentral.precompute.rna_type import insdc
from rnacentral_pipeline.rnacentral.precompute.rna_type import so_term


def rna_type_of(context, data) -> str:
    rna_type = so_term.rna_type_of(context, data)
    if not rna_type:
        insdc_rna_type = insdc.rna_type_of(data)
        if not insdc_rna_type:
            insdc_rna_type = 'ncRNA'
        rna_type = context.so_term_for(insdc_rna_type)
    return rna_type
