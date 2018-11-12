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

from . import families as fam

__doc__ = """
This module contains utility classes and functions for fetching and parsing
Rfam data from the FTP site.
"""


def name_to_insdc_type(handle):
    """
    Create a dict that maps from rfam family name (U1) to the INSDC RNA type
    (snRNA).
    """
    return {family.name: family.guess_insdc() for family in fam.parse(handle)}


def id_to_insdc_type(handle):
    """
    Create a dict that maps from rfam id (RF00003) to INSDC RNA type (snRNA).
    """
    return {family.id: family.guess_insdc() for family in fam.parse(handle)}


def id_to_pretty_name(handle):
    return {family.id: family.pretty_name for family in fam.parse(handle)}


def name_to_pretty_name(handle):
    return {family.name: family.pretty_name for family in fam.parse(handle)}


def name_to_suppression(handle):
    """
    Create a dict from the rfam family name (U1) to a flag indicating if the
    family data should be supressed in RNAcentral.
    """
    return {family.name: family.is_suppressed for family in fam.parse(handle)}


def id_to_suppression(handle):
    """
    Create a dict from the rfam family name (U1) to a flag indicating if the
    family data should be supressed in RNAcentral.
    """
    return {family.id: family.is_suppressed for family in fam.parse(handle)}
