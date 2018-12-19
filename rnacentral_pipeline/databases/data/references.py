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

import attr
from attr.validators import in_
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.helpers.hashes import md5

from . import utils

PMID_URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}+AND+SRC:MED&format=json'
DOI_URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{doi}+AND+SRC:MED&format=json'

KNOWN_SERVICES = {
    'doi',
    'pmid',
}


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

    authors = attr.ib(validator=is_a(basestring), converter=utils.optional_utf8)
    location = attr.ib(validator=is_a(basestring))
    title = attr.ib(
        validator=optional(is_a(basestring)),
        converter=utils.optional_utf8
    )
    pmid = attr.ib(validator=optional(is_a(int)))
    doi = attr.ib(validator=optional(is_a(basestring)))

    def md5(self):
        """
        Computes the MD5 hash of the reference.
        """
        return md5(''.join([
            (self.authors or ''),
            (self.location or ''),
            (self.title or ''),
        ]))

    def writeable_generic_pubmed(self):
        return [
            self.pmid,
            self.authors,
            self.location,
            self.title,
            self.doi,
        ]

    def writeable(self, accession):
        yield [
            self.md5(),
            accession,
            self.authors,
            self.location,
            self.title,
            self.pmid,
            self.doi,
        ]


@attr.s(frozen=True, hash=True)
class IdReference(object):
    namespace = attr.ib(validator=in_(KNOWN_SERVICES))
    external_id = attr.ib(validator=is_a(basestring))

    @classmethod
    def build(cls, ref_id):
        if isinstance(ref_id, int):
            return cls('pmid', str(ref_id))

        if isinstance(ref_id, basestring):
            if re.match('^\d+$', ref_id):
                return cls('pmid', ref_id)
            service, eid = ref_id.split(':', 1)
            service = service.lower()
            if service in KNOWN_SERVICES:
                return cls(service, eid)
        raise UnknownPublicationType(ref_id)

    @property
    def normalized_id(self):
        return '%s:%s' % (self.namespace, self.external_id)

    def external_url(self):
        if self.namespace == 'pmid':
            return PMID_URL.format(pmid=self.external_id)
        if self.namespace == 'doi':
            return DOI_URL.format(doi=self.external_id)
        raise ValueError("No URL for namespace %s" % self.namespace)

    def writeable(self, accession):
        yield [self.normalized_id, accession]

    def writeable_id(self):
        yield [self.normalized_id]
