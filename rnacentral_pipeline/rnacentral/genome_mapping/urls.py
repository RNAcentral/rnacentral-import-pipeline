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

import csv
import logging
import operator as op
import itertools as it
from ftplib import FTP
from contextlib import contextmanager


LOGGER = logging.getLogger(__name__)

try:
    import enum
except ImportError:
    from enum32 import enum


class NoTopLevelFiles(Exception):
    """
    Could not find any files that look like a genomic toplevel file to download.
    """
    pass


VERTEBRATES = {'ensembl', 'ensemblvertebrates', 'ensembl_vertebrates'}


@enum.unique
class FtpHost(enum.Enum):
    ensembl = 0
    ensembl_plants = 1
    ensembl_metazoa = 2
    ensembl_protists = 3
    ensembl_bacteria = 4
    ensembl_genomes = 5
    unknown = 6

    @classmethod
    def from_string(cls, raw):
        if raw.lower() in VERTEBRATES:
            return cls.ensembl
        if raw.lower() in {'ensemblgenomes'}:
            return cls.ensembl_genomes
        if hasattr(cls, raw.lower()):
            return getattr(cls, raw.lower())
        for name, member in cls.__members__.items():
            if name.replace('_', '') == raw.lower():
                return member
        raise ValueError("Unknown type of FTP host: %s" % raw)

    @property
    def host(self):
        if self is FtpHost.ensembl:
            return 'ftp.ensembl.org'
        if self is FtpHost.ensembl_plants:
            return 'ftp.ensemblgenomes.org'
        if self is FtpHost.ensembl_metazoa:
            return 'ftp.ensemblgenomes.org'
        if self is FtpHost.ensembl_protists:
            return 'ftp.ensemblgenomes.org'
        if self is FtpHost.ensembl_bacteria:
            return 'ftp.ensemblgenomes.org'
        if self is FtpHost.ensembl_genomes:
            return 'ftp.ensemblgenomes.org'
        raise ValueError("No host for %s" % self)

    def paths(self, species):
        if self is FtpHost.ensembl:
            return ['/pub/current_fasta/{species}/dna'.format(species=species)]
        else:
            template = '/pub/current/{group}/fasta/{species}/dna'
            if self is FtpHost.ensembl_plants:
                return [template.format(group='plants', species=species)]
            if self is FtpHost.ensembl_metazoa:
                return [template.format(group='metazoa', species=species)]
            if self is FtpHost.ensembl_protists:
                return [template.format(group='protists', species=species)]
            if self is FtpHost.ensembl_bacteria:
                return [template.format(group='bacteria', species=species)]
            if self is FtpHost.ensembl_genomes:
                gs = ['plants', 'metazoa', 'protists', 'bacteria']
                return [template.format(group=g, species=species) for g in gs]
        raise ValueError("No paths for %s" % self)


@contextmanager
def ftp(host):
    conn = FTP(host)
    conn.login()

    yield conn

    try:
        conn.quit()
    except Exception as err:
        LOGGER.info("Failed to close FTP connection")
        LOGGER.exception(err)


def toplevel_file(host, species, assembly_id, directory, files):
    upper_species = species[0].upper() + species[1:]
    base = '{species}.{assembly_id}.dna.{type}.fa.gz'.format(
        species=upper_species,
        assembly_id=assembly_id,
        type='{type}'
    )
    primary = base.format(type='primary_assembly')
    toplevel = base.format(type='toplevel')

    base_result = 'ftp://{host}{directory}/{file}'.format(
        directory=directory, 
        host=host,
        file='{file}'
    )
    if primary in files:
        return base_result.format(file=primary)
    elif toplevel in files:
        return base_result.format(file=toplevel)
    raise NoTopLevelFiles("%s %s" % (species, assembly_id))


def url_for(species, assembly_id, host=FtpHost.ensembl):

    if isinstance(host, str):
        return url_for(species, assembly_id, host=FtpHost.from_string(host))
    elif isinstance(host, FtpHost):
        if host is not FtpHost.unknown:
            with ftp(host.host) as conn:
                for path in host.paths(species):
                    try:
                        conn.cwd(path)
                        possible = set(conn.nlst())
                    except Exception:
                        continue
                    return toplevel_file(host.host, species, assembly_id, path, possible)

        else:
            for specific in list(FtpHost):
                if specific is FtpHost.unknown:
                    continue

                try:
                    return url_for(species, assembly_id, host=specific)
                except:
                    pass

        raise NoTopLevelFiles("%s %s" % (species, assembly_id))
    raise ValueError("Unknown type of Ftp Host: %s" % host)


def urls_for(handle):
    reader = csv.reader(handle)
    for (species, assembly_id, host) in reader:
        yield (species, assembly_id, url_for(species, assembly_id, host))


def write_urls_for(handle, output):
    writer = csv.writer(output)
    writer.writerows(urls_for(handle))
