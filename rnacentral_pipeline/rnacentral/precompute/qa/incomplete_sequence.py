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


FAMILIES_TO_CHECK = {
    'RF00001',  # 5S ribosomal RNA
    'RF00002',  # 5.8S ribosomal RNA
    'RF00005',  # tRNA
    'RF00177',  # Bacterial small subunit ribosomal RNA
    'RF01959',  # Archaeal small subunit ribosomal RNA
    'RF01960',  # Eukaryotic small subunit ribosomal RNA
    'RF02540',  # Archaeal large subunit ribosomal RNA
    'RF02541',  # Bacterial large subunit ribosomal RNA
    'RF02542',  # Microsporidia small subunit ribosomal RNA
    'RF02543',  # Eukaryotic large subunit ribosomal RNA
}


class Validator(object):
    name = 'incomplete_sequence'

    def status(self, _, data):
        """
        Detect if the given sequence is incomplete. This will work by checking if
        the Rfam hits cover enough of the sequence and the hit is complete enough.
        """

        hits = data.rfam_hits
        if len(hits) != 1:
            return False

        hit = hits[0]
        if hit.model not in FAMILIES_TO_CHECK:
            return False

        return hit.model_info.completeness <= 0.5 and \
            hit.sequence_info.completeness >= 0.9

    def message(self, _, data):
        name = data.rfam_hits[0].model_long_name
        url = data.rfam_hits[0].url
        message = 'Potential <a href="{url}">{name}</a> fragment'
        return message.format(name=name, url=url)
