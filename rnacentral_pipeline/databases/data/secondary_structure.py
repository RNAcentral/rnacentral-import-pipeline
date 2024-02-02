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
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.helpers.hashes import md5


@attr.s(frozen=True)
class SecondaryStructure(object):
    """
    This represents the secondary structure from GtRNAdb.
    """

    dot_bracket: str = attr.ib(validator=is_a(str))

    @classmethod
    def empty(cls):
        """
        Create an empty secondary structure.
        """
        return cls(dot_bracket="")

    def __bool__(self):
        """
        Check if this is empty.
        """
        return bool(self.dot_bracket)

    def __len__(self):
        return len(self.dot_bracket)

    @property
    def md5(self):
        """
        Compute the MD5 of the dot_bracket string.
        """
        return md5(self.dot_bracket)

    def writeable(self, accession):
        if not self.dot_bracket:
            return

        yield [
            accession,
            self.dot_bracket,
            self.md5,
        ]
