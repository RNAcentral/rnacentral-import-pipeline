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

import pytest

from tasks.config import ensembl

from tasks.ensembl.utils import generic


@pytest.mark.parametrize('species,state', [
    ('Mus_musculus_akrj', False),
    ('mus_musculus_akrj', False),
    ('saccharomyces_cerevisiae', False),
    ('caenorhabditis_elegans', False),
    ('drosophila_melanogaster', False),
    ("dipodomys_ordii", True),
    ("echinops_telfairi", True),
    ("equus_caballus", True),
    ("erinaceus_europaeus", True),
    ("felis_catus", True),
    ("ficedula_albicollis", True),
    ("fukomys_damarensis", True),
    ("gadus_morhua", True),
    ("gallus_gallus", True),
    ("gasterosteus_aculeatus", True),
    ("gorilla_gorilla", True),
    ("heterocephalus_glaber_female", True),
    ("heterocephalus_glaber_male", True),
    ("homo_sapiens", True),
    ("ictidomys_tridecemlineatus", True),
    ("jaculus_jaculus", True),
    ("latimeria_chalumnae", True),
    ("lepisosteus_oculatus", True),
    ("loxodonta_africana", True),
    ("macaca_mulatta", True),
    ("meleagris_gallopavo", True),
    ("mesocricetus_auratus", True),
    ("microcebus_murinus", True),
    ("microtus_ochrogaster", True),
    ("monodelphis_domestica", True),
    ("mus_caroli", True),
    ("mus_musculus", True),
    ("mus_musculus_129s1svimj", False),
    ("mus_musculus_aj", False),
    ("mus_musculus_akrj", False),
    ("mus_musculus_balbcj", False),
    ("mus_musculus_c3hhej", False),
    ("mus_musculus_c57bl6nj", False),
    ("mus_musculus_casteij", False),
    ("mus_musculus_dba2j", False),
    ("mus_musculus_fvbnj", False),
    ("mus_musculus_lpj", False),
    ("mus_musculus_nodshiltj", False),
    ("mus_musculus_nzohlltj", False),
    ("mus_musculus_pwkphj", False),
    ("mus_musculus_wsbeij", False),
    ("mus_pahari", True),
    ("mus_spretus_spreteij", True),
    ("mustela_putorius_furo", True),
    ("myotis_lucifugus", True),
    ("nannospalax_galili", True),
    ("nomascus_leucogenys", True),
    ("notamacropus_eugenii", True),
    ("ochotona_princeps", True),
    ("octodon_degus", True),
    ("oreochromis_niloticus", True),
    ("ornithorhynchus_anatinus", True),
])
def test_it_knows_which_species_are_included_by_default(species, state):
    config = ensembl()
    assert generic.allowed_species(config, species) == state


@pytest.mark.parametrize('species,state', [
    ('Mus_musculus_akrj', True),
    ('mus_musculus_akrj', True),
    ("dipodomys_ordii", False),
    ("echinops_telfairi", False),
    ("equus_caballus", False),
    ("mus_musculus", False),
    ("mus_pahari", False),
    ("mus_spretus_spreteij", False),
])
def test_it_knows_if_something_is_mouse_strain(species, state):
    assert generic.is_mouse_strain(species) == state


@pytest.mark.parametrize('species', [
    ("loxodonta_africana"),
    ("macaca_mulatta"),
    ("meleagris_gallopavo"),
    ("mesocricetus_auratus"),
    ("microcebus_murinus"),
    ("microtus_ochrogaster"),
    ("monodelphis_domestica"),
    ("mus_caroli"),
    ("mus_musculus"),
    ("mus_musculus_129s1svimj"),
    ("mus_musculus_aj"),
    ("mus_musculus_akrj"),
    ("mus_musculus_balbcj"),
    ("mus_musculus_c3hhej"),
    ("mus_musculus_c57bl6nj"),
    ("mus_musculus_casteij"),
    ("mus_musculus_dba2j"),
    ("mus_musculus_fvbnj"),
    ("mus_musculus_lpj"),
    ("mus_musculus_nodshiltj"),
    ("mus_musculus_nzohlltj"),
    ("mus_musculus_pwkphj"),
    ("mus_musculus_wsbeij"),
])
def test_can_be_configured_to_allow_mouse_strains(species):
    config = ensembl()
    config.allow_model_organisms = True
    config.exclude_mouse_strains = False
    assert generic.allowed_species(config, species) is True


@pytest.mark.parametrize('species,state', [
    ('Mus_musculus_akrj', False),
    ('mus_musculus_akrj', False),
    ("mus_musculus_casteij", False),
    ('saccharomyces_cerevisiae', True),
    ('caenorhabditis_elegans', True),
    ('drosophila_melanogaster', True),
    ("mus_spretus_spreteij", True),
    ("mustela_putorius_furo", True),
    ("myotis_lucifugus", True),
])
def test_alloweed_species_can_be_configured_to_allow_models(species, state):
    config = ensembl()
    config.allow_model_organisms = True
    config.exclude_mouse_strains = True
    assert generic.allowed_species(config, species) == state


@pytest.mark.parametrize('species', [
    ('Mus_musculus_akrj'),
    ('mus_musculus_akrj'),
    ("mus_musculus_casteij"),
    ('saccharomyces_cerevisiae'),
    ('caenorhabditis_elegans'),
    ('drosophila_melanogaster'),
    ("mus_spretus_spreteij"),
    ("mustela_putorius_furo"),
    ("myotis_lucifugus"),
])
def test_allowed_species_can_be_configured_to_both(species):
    config = ensembl()
    config.allow_model_organisms = True
    config.exclude_mouse_strains = False
    assert generic.allowed_species(config, species) is True
