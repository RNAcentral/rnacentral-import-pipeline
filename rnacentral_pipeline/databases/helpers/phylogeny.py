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

import logging
import typing as ty
from functools import lru_cache
from time import sleep

import requests
import simplejson

TAX_URL = "https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{taxon_id}"

SPECIES_URL = "https://www.ebi.ac.uk/ena/taxonomy/rest/any-name/{species}"

LOGGER = logging.getLogger(__name__)


class UnknownTaxonId(Exception):
    """
    This is raised if we cannot get the required information about the given
    taxon id. This happens with the taxon id is obsoleted, or invalid.
    """

    pass


class FailedTaxonId(Exception):
    """
    This is raised when we cannot fetch information on a taxon id.
    """

    pass


@lru_cache()
def phylogeny(taxon_id: int) -> ty.Dict[str, str]:
    """
    Call the EBI taxonomy API to get the phylogenetic information for the given
    taxon id. This will cache requests to the same taxon id. This will retry
    the request up to 10 times (with a sleep between each one) in the case of
    500 errors. These seem to happen sometimes but go away with more requests.
        However in the case of 400 errors this will fail on the first attempt.
    """

    for count in range(10):
        response = requests.get(TAX_URL.format(taxon_id=taxon_id))
        try:
            response.raise_for_status()
            data = response.json()
            break
        except simplejson.errors.JSONDecodeError:
            sleep(0.15 * (count + 1) ** 2)
            continue
        except requests.HTTPError as err:
            if response.status_code == 500:
                sleep(0.15 * (count + 1) ** 2)
                continue
            elif response.status_code == 404:
                raise UnknownTaxonId(taxon_id)
            else:
                LOGGER.exception(err)
                raise FailedTaxonId("Unknown error")
    else:
        raise FailedTaxonId("Could not get taxon id for %s" % taxon_id)

    if not data:
        raise FailedTaxonId("Somehow got no data")

    return data


def lineage(taxon_id: int) -> str:
    """
    Call the EBI taxonomy API to fetch the lineage for the given
    """

    data = phylogeny(taxon_id)
    return str(
        "{lineage}{name}".format(
            lineage=data.get("lineage", None), name=data["scientificName"]
        )
    )


def common_name(taxon_id: int) -> ty.Optional[str]:
    """
    Get the common name, if any for the given taxon id. If no common name
    exists then None is returned.
    """
    return phylogeny(taxon_id).get("commonName", None)


def species(taxon_id: int) -> str:
    """
    Get a standardized species name for the given taxon id.
    """

    data = phylogeny(taxon_id)
    return data["scientificName"]


def division(taxon_id: int) -> str:
    """
    Get the annotated division for the given taxon id.
    """

    data = phylogeny(taxon_id)
    return data["division"]


@lru_cache
def taxid(species: str) -> int:
    """
    Get the taxid for a given species
    Re-use request logic from phylogeny, but this uses a different endpoint, so
    can't directly reuse
    """

    for count in range(10):
        response = requests.get(SPECIES_URL.format(species=species))
        try:
            response.raise_for_status()
            data = response.json()
            break
        except simplejson.errors.JSONDecodeError:
            sleep(0.15 * (count + 1) ** 2)
            continue
        except requests.HTTPError as err:
            if response.status_code == 500:
                sleep(0.15 * (count + 1) ** 2)
                continue
            elif response.status_code == 404:
                raise UnknownTaxonId(taxon_id)
            else:
                LOGGER.exception(err)
                raise FailedTaxonId("Unknown error")
    else:
        raise FailedTaxonId("Could not get taxon id for %s" % species)

    if not data:
        raise FailedTaxonId("Somehow got no data")
    return int(data[0]["taxId"])
