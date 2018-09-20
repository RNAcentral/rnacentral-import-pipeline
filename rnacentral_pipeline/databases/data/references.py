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

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.helpers.hashes import md5

from . import utils


@attr.s(frozen=True)
class Reference(object):
    """
    This stores the data for a reference that will be written to out to csv
    files.
    """

    authors = attr.ib(validator=is_a(basestring), convert=utils.optional_utf8)
    location = attr.ib(validator=is_a(basestring))
    title = attr.ib(
        validator=optional(is_a(basestring)),
        convert=utils.optional_utf8
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
