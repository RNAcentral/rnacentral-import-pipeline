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

import re
import enum

import attr
from attr.validators import in_
from attr.validators import optional
from attr.validators import instance_of as is_a

from furl import furl

from rnacentral_pipeline.databases.helpers.hashes import md5

from . import utils


@enum.unique
class KnownServices(enum.Enum):
    doi = 0
    pmid = 1
    pmcid = 2

    @classmethod
    def from_name(cls, name):
        return getattr(cls, name.lower())


class UnknownPublicationType(Exception):
    """
    Raised when we are trying to build a reference but it is to an unknown type
    of identifier.
    """
    pass


@attr.s(frozen=True)
class Reference(object):
    """
    This stores the data for a reference that will be written to out to csv
    files.
    """

    authors = attr.ib(validator=is_a(str), converter=str)
    location = attr.ib(validator=is_a(str), converter=str)
    title = attr.ib(validator=optional(is_a(str)), converter=str)
    pmid = attr.ib(validator=optional(is_a(int)))
    doi = attr.ib(validator=optional(is_a(str)))
    pmcid = attr.ib(validator=optional(is_a(str)), default=None)

    def md5(self):
        """
        Computes the MD5 hash of the reference.
        """

        title = self.title if self.title else ''
        data = [
            self.authors,
            self.location,
            title,
        ]

        return md5(''.join(data).encode('utf-8'))

    def writeable_generic_pubmed(self):
        return [
            self.pmid,
            self.authors,
            self.location,
            self.title,
            self.doi,
        ]

    def writeable(self, extra):
        data = [self.md5()]
        if isinstance(extra, str):
            data.append(extra)
        else:
            data.extend(extra)

        rest = [
            self.authors,
            self.location,
            self.title,
            self.pmid,
            self.doi,
        ]
        yield data + rest

    @property
    def id_reference(self):
        if self.pmid:
            return IdReference(
                namespace=KnownServices.pmid, 
                external_id=str(self.pmid),
            )
        if self.doi:
            return IdReference(
                namespace=KnownServices.doi, 
                external_id=self.doi,
            )
        if self.pmcid:
            return IdReference(
                namespace=KnownServices.pmcid, 
                external_id=self.pmcid,
            )
        raise ValueError("Cannot build IdReference for %s" % self)

    @property
    def id_references(self):
        refs = set()
        for namespace in KnownServices:
            value = getattr(self, namespace.name)
            if not value:
                continue
            eid = str(value)
            refs.add(IdReference(namespace=namespace, external_id=str(eid)))
        return refs


@attr.s(frozen=True, hash=True)
class IdReference(object):
    namespace = attr.ib(validator=is_a(KnownServices))
    external_id = attr.ib(validator=is_a(str)) 

    @classmethod
    def build(cls, ref_id):
        if isinstance(ref_id, int):
            return cls(KnownServices.pmid, str(ref_id))

        if not isinstance(ref_id, str):
            raise UnknownPublicationType(ref_id)

        ref_id = str(ref_id.strip())
        if re.match(r'^\d+$', ref_id):
            return cls(KnownServices.pmid, ref_id)

        if re.match(r'^PMC\d+$', ref_id, re.IGNORECASE):
            return cls(KnownServices.pmcid, ref_id.upper())

        if ':' not in ref_id:
            raise UnknownPublicationType("Could not parse: " + ref_id)

        service, eid = ref_id.split(':', 1)
        service = service.lower()
        service = KnownServices.from_name(service)
        if service is KnownServices.pmcid:
            eid = eid.upper()
            if not eid.startswith('PMC'):
                eid = 'PMC' + eid
        return cls(service, eid)

    @property
    def normalized_id(self):
        return '%s:%s' % (self.namespace.name, self.external_id)

    def external_url(self):
        base = furl('https://www.ebi.ac.uk/europepmc/webservices/rest/search')
        base.args['format'] = 'json'
        base.args['pageSize'] = 1000

        if self.namespace is KnownServices.pmid:
            query = '{pmid} AND SRC:MED'.format(pmid=self.external_id)
            base.args['query'] = query
            return base.url

        if self.namespace is KnownServices.doi or \
                self.namespace is KnownServices.pmcid:
            suffix = self.external_id
            if self.namespace is KnownServices.doi:
                suffix = '"%s"' % suffix
            base.args['query'] = '%s:%s' % (self.namespace.name.upper(), suffix)
            return base.url

        raise ValueError("No URL for namespace %s" % self.namespace)

    def writeable(self, accession):
        yield [self.normalized_id, accession]

    def writeable_id(self):
        yield [self.normalized_id]
