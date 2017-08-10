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
import logging

import attr
from attr.validators import instance_of as is_a

from tasks.utils.ftp import ftp_context

LOGGER = logging.getLogger(__name__)


@attr.s(frozen=True)  # pylint: disable=R0903
class SpeciesDescription(object):
    """
    A simple class to store the information about the species data to import
    from Ensembl.

    Properties
    ----------
    species_name : str
        The name of the species like 'Homo sapiens'.

    filenames : list
        List of filenames (not full paths) of EMBL file to import.
    """
    species_name = attr.ib(validator=is_a(basestring))
    filenames = attr.ib(validator=is_a(list))


def is_data_file(name):
    """
    Check if the given filename contains genome data to import. By default
    this will only use data from chromosomes. However there is an option to
    allow data from non-chromosomal files. In Ensembl's FTP site each
    species contains utility files like 'README' as well as EMBL files that
    contain data to import. This method helps distinguish the two.

    Parameters
    ----------
    name : str
    The filename, may be a full path as well.

    allow_nonchromosomal : bool, False
    If this should accept EMBL files which are for non-chromosomal data
    as well as chromosomal data.

    Returns
    -------
    is_data_file : bool
    True if the filename is for data this can import.
    """

    filename = os.path.basename(name)
    return 'chromosome' in filename or 'nonchromosomal' in filename


def release_to_path(release):
    """
    Turn an ensembl release id to the path on the FTP site for data from that
    release. The release can be a number like 99, or the string 'current'.

    Returns
    -------
    path : str
        The path.
    """

    if release == 'current':
        return 'pub/current_embl'
    return 'pub/release-{release}/embl/'.format(release=release)


def allowed_species(config, name):
    """
    Check if given the current configuration if the given species name is one
    that should be imported.
    """

    if config.allow_model_organisms:
        return True
    model_organisms = config.model_organism_set()
    return name.lower().replace(' ', '_') not in model_organisms


def description_of(ftp, name):
    """
    Compute the SpeciesDescription for the species with the given name.
    This will build a description using only chromosomal files, unless
    there are no such files in which case the non-chromosomal files will be
    used. If there are still no files found then None is returned.

    Parameters
    ----------
    name : str
        Name of the species in Ensembl to import.

    Returns
    -------
    description : SpeciesDescription or None
        The description, if any of the given species.
    """

    cleaned = os.path.basename(name).lower().replace('_', ' ')
    cleaned = cleaned[0].upper() + cleaned[1:]
    path = name.lower()
    names = [f for f in ftp.nlst(path) if is_data_file(f)]
    if names:
        return SpeciesDescription(
            species_name=cleaned,
            filenames=names,
        )
    LOGGER.info("No importable data found for %s", name)
    return None


def known_species(config, names):
    """
    Get the specified descriptions. If the given pattern is a name the
    species with that name will be selected.

    Parameters
    ----------
    names : list
        Names of the species to get a description for. If given 'all' then all
        known Ensembl species will be used.

    Returns:
    --------
    descriptions : list
        A list of either a SpeciesDescription object or None.
    """

    host = config.ftp_host
    base = release_to_path(config.release)
    with ftp_context(host, base=base) as ftp:
        files = []
        if names == 'all':
            files = ftp.nlst()
        else:
            files = sorted(names)

        descriptions = []
        for name in files:
            if not allowed_species(config, name):
                continue
            description = description_of(ftp, name)
            if description:
                descriptions.append(description)

        if not description:
            raise ValueError("No suitable organisms found")
        return descriptions
