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

EXPECTED_MATCHES = {
    'rRNA': {
        'RF00001',  # 5S ribosomal RNA
        'RF00002',  # 5.8S ribosomal RNA
        'RF00177',  # Bacterial small subunit ribosomal RNA
        'RF01959',  # Archaeal small subunit ribosomal RNA
        'RF01960',  # Eukaryotic small subunit ribosomal RNA
        'RF02540',  # Archaeal large subunit ribosomal RNA
        'RF02541',  # Bacterial large subunit ribosomal RNA
        'RF02542',  # Microsporidia small subunit ribosomal RNA
        'RF02543',  # Eukaryotic large subunit ribosomal RNA
        'RF02547',  # mito 5S RNA
    },
    'tRNA': {
        'RF00005',  # tRNA
        'RF01852',  # Selenocysteine tRNA
    },
}


class Validator(object):
    name = 'missing_rfam_match'

    def status(self, rna_type, data):
        """
        Detect if the given sequence, with the given RNA type is missing an Rfam
        match.
        """

        if rna_type not in EXPECTED_MATCHES:
            return False

        if rna_type == 'tRNA' and \
                data.is_mitochondrial() and \
                not data.rfam_hits:
            return False

        required = EXPECTED_MATCHES[rna_type]
        hits = {h.model for h in data.rfam_hits}
        return not hits.intersection(required)

    def href(self, model_id):
        link = '<a href="http://rfam.org/family/{model_id}">{model_id}</a>'
        return link.format(model_id=model_id)

    def message(self, rna_type, _):
        possible = sorted(EXPECTED_MATCHES[rna_type])

        article = 'the'
        if len(possible) > 1:
            article = 'a'

        models = [self.href(p) for p in sorted(possible)]
        raw = 'No match to {article} {rna_type} Rfam model ({possible})'
        return raw.format(
            rna_type=rna_type,
            article=article,
            possible=', '.join(models),
        )
