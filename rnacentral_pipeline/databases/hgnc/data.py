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

from rnacentral_pipeline.databases.data import Entry


@attr.s(slots=True, frozen=True)
class MappableEntry(object):
    attribute_name = attr.ib(validator=is_a(basestring))
    value = attr.ib()
    sequence = attr.ib(validator=is_a(basestring))


@attr.s()
class KnownMapper(object):
    entries = attr.ib()

    def matching_updates(self, partial):
        pass

    def store_all(self, given):
        pass

    def as_entries(self, partial):
        """
        Turn a PartialEntry into a full Entry by finding the matching update
        and apply it.
        """

        data = []
        for update in self.matching_updates(partial):
            extra_fields = set(f.name for f in attr.fields(update.__class__))
            for field in attr.fields(Entry):
                target = partial
                if field.name in extra_fields:
                    target = update
                data.append(getattr(target, field.name))
            yield Entry(*data)  # pylint: disable=star-args
