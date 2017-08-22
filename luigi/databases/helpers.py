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

import hashlib
from time import sleep
import logging

import requests

from functools32 import lru_cache

TAX_URL = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{taxon_id}'

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
def phylogney(taxon_id):
    """
    Call the EBI taxonomy API to get the phylogenetic information for the given
    taxon id. This will cache requests to the same taxon id. This will retry
    the request up to 10 times (with a sleep between each one) in the case of
    500 errors. These seem to happen sometimes but go away with more requests.
        However in the case of 400 errors this will fail on the first attempt.
    """

    for count in xrange(10):
        response = requests.get(TAX_URL.format(taxon_id=taxon_id))
        try:
            response.raise_for_status()
            break
        except requests.HTTPError as err:
            if response.status_code == 500:
                sleep(0.1 * (count + 1))
                continue
            else:
                print(err)
                raise UnknownTaxonId(taxon_id)
    else:
        raise FailedTaxonId("Could not get taxon id for %s" % taxon_id)

    data = response.json()
    assert data, "Somehow got no data"
    return data


def lineage(taxon_id):
    """
    Call the EBI taxonomy API to fetch the lineage for the given
    """

    data = phylogney(taxon_id)
    return '{lineage}{name}'.format(
        lineage=data['lineage'],
        name=data['scientificName']
    )


def common_name(taxon_id):
    """
    Get the common name, if any for the given taxon id. If no common name
    exists then None is returned.
    """

    data = phylogney(taxon_id)
    return data.get('common_name', None)


def species(taxon_id):
    """
    Get a standardized species name for the given taxon id.
    """

    data = phylogney(taxon_id)
    return data['scientificName']


def current_taxon_id(taxon_id):
    """
    If the taxon id is oudated then this will look up the current taxon id,
    otherwise it will return the given taxon id.
    """
    return taxon_id


def md5(data):
    """
    Get the MD5 hash as a string.
    """
    return hashlib.md5(data).hexdigest()


def crc64(input_string):
    """
    Python re-implementation of SWISS::CRC64
    Adapted from:
    http://code.activestate.com/recipes/259177-crc64-calculate-the-cyclic-redundancy-check/
    """
    POLY64REVh = 0xd8000000L
    CRCTableh = [0] * 256
    CRCTablel = [0] * 256
    isInitialized = False
    crcl = 0
    crch = 0
    if isInitialized is not True:
        isInitialized = True
        for i in xrange(256):
            partl = i
            parth = 0L
            for _ in xrange(8):
                rflag = partl & 1L
                partl >>= 1L
                if parth & 1:
                    partl |= (1L << 31L)
                parth >>= 1L
                if rflag:
                    parth ^= POLY64REVh
            CRCTableh[i] = parth
            CRCTablel[i] = partl

    for item in input_string:
        shr = 0L
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8L
        temp1l = (crcl >> 8L) | shr
        tableindex = (crcl ^ ord(item)) & 0xFF

        crch = temp1h ^ CRCTableh[tableindex]
        crcl = temp1l ^ CRCTablel[tableindex]
    return "%08X%08X" % (crch, crcl)
