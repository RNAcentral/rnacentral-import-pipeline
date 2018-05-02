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


EXCLUDED_TERMS = {
    'GO:0008049',
    'GO:0042981',
    'GO:0042749',
    'GO:0050789',
    'GO:0006810',
}

EXCLUDED_MIRNA = {'GO:0035068', 'GO:0006396'}


GO_REPLACEMENTS = {
    'GO:0044005': 'GO:0051819',
}


def go_term_mapping(version='CURRENT'):
    for family in load_families(version=version):
        for (go_term_id, _) in family.go_terms:
            if go_term_id in EXCLUDED_TERMS:
                continue

            if family.rna_type == 'Gene; miRNA' and go_term_id in EXCLUDED_MIRNA:
                continue

            go_term_id = GO_REPLACEMENTS.get(go_term_id, go_term_id)

            yield {
                'go_term_id': go_term_id,
                'rfam_model_id': family.id,
            }
