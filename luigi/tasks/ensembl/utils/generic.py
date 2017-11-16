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

import os
import re

from databases.ensembl.parsers import EnsemblParser
from databases.ensembl.parsers import GencodeParser


def is_gencode_file(config, input_file):
    """
    Check if the given file contains GENCODE data. This is done by checking if
    file starts with a species name that is in the configured gencode species.
    """

    species = os.path.basename(input_file).split('.')[0]
    normalized = normalize_species_name(species)
    return normalized in config.gencode_species_set()


def parser_class(config, input_file):
    """
    Given a configuration object and the input filename this will determine
    which parser class to use. If it is a GENCODE organism then the GENCODE
    paresr will be returned, otherwise the more generic one will be provided.
    """

    if is_gencode_file(config, input_file):
        return GencodeParser()
    return EnsemblParser()


def normalize_species_name(species):
    """
    This will put species names into a standard format. That is lower case,
    without leading or trailing whitespace and with spaces replaced by '_'.
    """
    return species.strip().lower().replace(' ', '_')


def is_mouse_strain(species):
    """
    Check if a species name is likely a mouse strain.
    """
    return bool(re.search('^mus_musculus_.*', species, re.IGNORECASE))


def allowed_species(config, species):
    """
    Check if given the current configuration if the given species name is one
    that should be imported.
    """

    if config.allow_model_organisms:
        return not (config.exclude_mouse_strains and is_mouse_strain(species))

    if normalize_species_name(species) in config.model_organism_set():
        return False
    return not (config.exclude_mouse_strains and is_mouse_strain(species))
