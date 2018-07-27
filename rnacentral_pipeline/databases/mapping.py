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

from collections import defaultdict
from collections import OrderedDict

import attr
from attr.validators import instance_of as is_a

from .data import Entry


def class_without_fields(cls, *names, **kwargs):
    """
    Build a new class that has all the same fields as a given class, excluding
    the listed fields. The given keywords will be used as options when building
    the class.
    """

    fields = OrderedDict()
    ignore = set(names)
    for field in attr.fields(cls):
        if field.name in ignore:
            continue

        fields[field.name] = attr.ib(
            default=field.default,
            validator=field.validator,
            repr=field.repr,
            cmp=field.cmp,
            hash=field.hash,
            init=field.init,
            convert=field.convert,
            metadata={'name': field.name},
        )

    return attr.make_class('Simplified' + cls.__name__, fields, **kwargs)


SequencelessEntries = class_without_fields(Entry, 'sequence', frozen=True)


@attr.s(slots=True, frozen=True)
class EntryCompleter(object):
    matchable_id = attr.ib(validator=is_a(basestring))
    delta = attr.ib(validator=is_a(dict))

    @classmethod
    def from_dict(cls, raw):
        matchable_id = raw.pop('id')
        return cls(
            matchable_id=matchable_id,
            delta=dict(raw),
        )

    def as_entry(self, partial):
        data = []
        for field in attr.fields(Entry):
            if field.name in self.delta:
                data.append(self.delta[field.name])
            else:
                data.append(getattr(partial, field.name))
        return Entry(*data)  # pylint: disable=star-args


@attr.s(frozen=True)
class Matcher(object):
    database_name = attr.ib(validator=is_a(basestring))
    updates = attr.ib(validator=is_a(dict))

    @classmethod
    def from_iterable(cls, key, iterable):
        return cls(
            database_name=key,
            updates=[EntryCompleter.from_dict(val) for val in iterable],
        )

    def store(self, completer):
        self.updates[completer.value].append(completer)

    def matching_updates(self, partial):
        value = partial.xref_data.get(self.database_name)
        return self.updates.get(value, [])

    def as_entries(self, partial):
        for update in self.matching_updates(partial):
            yield update.as_entry(partial)


@attr.s()
class MultiMatcher(object):
    matchers = attr.ib()

    @classmethod
    def build(cls, *matchers):
        return cls(list(matchers))

    def as_entries(self, partial):
        for matcher in self.matchers:
            if matcher.matching_updates(partial):
                for entry in matcher.as_entries(partial):
                    yield entry
                break
