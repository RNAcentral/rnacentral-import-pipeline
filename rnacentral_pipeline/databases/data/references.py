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

import six

import attr
from attr.validators import in_
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.helpers.hashes import md5

from . import utils

PMID_URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}+AND+SRC:MED&format=json'
DOI_URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{doi}+AND+SRC:MED&format=json'
PMC_URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={pmc}+AND+SRC:MED&format=json'

KNOWN_SERVICES = {
    'doi',
    'pmid',
    'pmcid',
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

    authors = attr.ib(validator=is_a(six.text_type))
    location = attr.ib(validator=is_a(six.text_type))
    title = attr.ib(validator=optional(is_a(six.text_type)))
    pmid = attr.ib(validator=optional(is_a(int)))
    doi = attr.ib(validator=optional(is_a(six.text_type)))
    pmcid = attr.ib(validator=optional(is_a(six.text_type)), default=None)

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
        if isinstance(extra, six.string_types):
            data.append(extra)
        else:
            data.extend(extra)

        data.extend([
            self.authors.encode('utf-8'),
            self.location.encode('utf-8'),
            self.title.encode('utf-8'),
            self.pmid,
            self.doi.encode('utf-8'),
        ])

        yield data


@attr.s(frozen=True, hash=True)
class IdReference(object):
    namespace = attr.ib(validator=in_(KNOWN_SERVICES))
    external_id = attr.ib(validator=is_a(six.text_type))

    @classmethod
    def build(cls, ref_id):
        if isinstance(ref_id, int):
            return cls('pmid', six.text_type(ref_id))

        if isinstance(ref_id, six.text_type):
            ref_id = six.text_type(ref_id.strip())
            if re.match(r'^\d+$', ref_id):
                return cls('pmid', ref_id)
            if re.match(r'^PMC\d+', ref_id, re.IGNORECASE):
                return cls('pmcid', ref_id.upper())
            if ':' not in ref_id:
                raise UnknownPublicationType("Could not parse: " + ref_id)
            service, eid = ref_id.split(':', 1)
            service = service.lower()
            if service in KNOWN_SERVICES:
                if service == 'pmcid':
                    eid = eid.upper()
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
        if self.namespace == 'pmcid':
            return PMC_URL.format(pmc=self.external_id)
        raise ValueError("No URL for namespace %s" % self.namespace)

    def writeable(self, accession):
        yield [self.normalized_id, accession]

    def writeable_id(self):
        yield [self.normalized_id]
