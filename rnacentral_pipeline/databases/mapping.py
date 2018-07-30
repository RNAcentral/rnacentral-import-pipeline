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


@attr.s(frozen=True)
class SequencelessEntry(class_without_fields(Entry, 'sequence')):
    """
    This represents an Entry object which has all possible fields, execpt the
    'sequence' field. Generally when mapping sequences we know a lot of
    information about the Entry execpt the Sequence. This provides a simple way
    of representing such cases. To ensure that these data do not get written by
    accident the is_valid method will always return False, and the writing
    methods will all raise an Exception.
    """

    def is_valid(self):
        return False

    def write_ac_info(self):
        raise ValueError("Not possible for SequencelessEntries")

    def write_secondary_structure(self):
        raise ValueError("Not possible for SequencelessEntries")

    def write_seq_long(self):
        raise ValueError("Not possible for SequencelessEntries")

    def write_seq_short(self):
        raise ValueError("Not possible for SequencelessEntries")

    def write_refs(self):
        raise ValueError("Not possible for SequencelessEntries")

    def write_genomic_locations(self):
        raise ValueError("Not possible for SequencelessEntries")


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
        updates = defaultdict(list)
        for val in iterable:
            update = EntryCompleter.from_dict(val)
            updates[update.matchable_id].append(update)

        return cls(
            database_name=key,
            updates=dict(updates),
        )

    def is_matchable(self, partial):
        keys = partial.xref_data.get(self.database_name, [])
        return any(k in self.updates for k in keys)

    def store(self, completer):
        self.updates[completer.value].append(completer)

    def matching_updates(self, partial):
        keys = partial.xref_data.get(self.database_name, [])
        updates = []
        for key in keys:
            # pylint: disable=no-member
            updates.extend(self.updates.get(key, []))
        return updates

    def as_entries(self, partial):
        for update in self.matching_updates(partial):
            yield update.as_entry(partial)


@attr.s()
class MultiMatcher(object):
    matchers = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, *matchers):
        return cls(list(matchers))

    def is_matchable(self, partial):
        return any(m.is_matchable(partial) for m in self.matchers)

    def matching_database_name(self, partial):
        for matcher in self.matchers:
            if matcher.is_matchable(partial):
                return matcher.database_name
        return None

    def as_entries(self, partial):
        for matcher in self.matchers:
            if matcher.matching_updates(partial):
                for entry in matcher.as_entries(partial):
                    yield entry
                break
