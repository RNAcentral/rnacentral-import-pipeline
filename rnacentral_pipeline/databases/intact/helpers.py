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

from rnacentral_pipeline.databases import data


def interactor_from_rnacentral(interactor):
    possible = [interactor.id] + interactor.alt_ids + interactor.aliases
    return any(i.key == 'rnacentral' for i in possible)


def involves_rnacentral(interaction: data.Interaction):
    return interactor_from_rnacentral(interaction.interactor1) or \
        interactor_from_rnacentral(interaction.interactor2)
