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

from tasks.ensembl.utils.generic import allowed_species
from tasks.ensembl.utils.generic import normalize_species_name
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


def data_files_path(release, datatype):
    """
    Get the path on the remote server to the given type of data files for the
    given release. This does not check if such a path exists, only constructs
    what the path will be. Release should be either 'current' or a known
    release number.
    """

    if release == 'current':
        return 'pub/current_%s' % datatype
    return 'pub/release-{release}/{datatype}/'.format(
        release=release,
        datatype=datatype,
    )


def embl_directory(release):
    """
    Turn an ensembl release id to the path on the FTP site for data from that
    release. The release can be a number like 99, or the string 'current'.

    Returns
    -------
    path : str
        The path.
    """
    return data_files_path(release, 'embl')


def species_file_path(config, filename, datatype='embl'):
    """
    Given a filename return the path to where that file will be in the Ensembl
    FTP site.
    """

    return '{base}/{filename}'.format(
        base=data_files_path(config.release, datatype),
        filename=filename,
    )


def known_species(config):
    """
    Get the specified descriptions. If the given pattern is a name the
    species with that name will be selected.

    Returns:
    --------
    descriptions : list
        A list of either a SpeciesDescription object or None.
    """

    host = config.ftp_host
    base = embl_directory(config.release)
    configured = config.species_names()
    with ftp_context(host, base=base) as ftp:
        files = []
        if configured == 'all':
            files = ftp.nlst()
        else:
            files = sorted(configured)

        names = []
        for name in files:
            if not allowed_species(config, name):
                continue
            path = normalize_species_name(name)
            if any(f for f in ftp.nlst(path) if is_data_file(f)):
                names.append(name)

    if not names:
        raise ValueError("No organisms found for: %s" % configured)
    return names


def species_description(config, species):
    """
    Generate a description of the species. This will be a SpeciesDescription
    object.
    """

    base = embl_directory(config.release)
    with ftp_context(config.ftp_host, base=base) as ftp:
        path = normalize_species_name(species)
        names = [f for f in ftp.nlst(path) if is_data_file(f)]

    if names:
        return SpeciesDescription(species_name=species, filenames=names)

    LOGGER.info("No importable data found for %s", species)
    return None
