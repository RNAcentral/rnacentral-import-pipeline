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

from rnacentral_pipeline import psql


@attr.s()
class QaStatus(object):
    """
    This represents an update to the QA status table.
    """

    incomplete_sequence = attr.ib(validator=is_a(bool))
    possible_contamination = attr.ib(validator=is_a(bool))
    missing_rfam_match = attr.ib(validator=is_a(bool))
    is_repetitive = attr.ib(validator=is_a(bool))
    mismatching_rna_type = attr.ib(validator=is_a(bool))
    messages = attr.ib(validator=is_a(list), default=attr.Factory(list))
    no_data = attr.ib(validator=optional(is_a(bool)), default=False)

    @classmethod
    def from_validators(cls, validators, *args, **kwargs):
        """
        This will build a new QaStatus object update given the particular RNA
        type and a Sequence object.
        """
        status = {'messages': [], 'no_data': False}
        for validator in validators:
            current = validator.status(*args, **kwargs)
            status[validator.name] = current
            if current:
                status['messages'].append(validator.message(*args, **kwargs))
        return cls(**status)

    @classmethod
    def empty(cls):
        """
        Build a QaStatus that will not produce any writeable data. Objects
        produced by this should not be used for any updates.
        """
        return cls(False, False, False, False, False, no_data=True)

    @property
    def has_issue(self):
        """
        Check if this QA update indicates if there is any issue.
        """

        return (
            self.incomplete_sequence or
            self.possible_contamination or
            self.missing_rfam_match
        )

    def writeable(self, upi, taxid):
        """
        Create a writeable array for writing CSV files.
        """

        if self.no_data:
            return None

        return [
            '%s_%i' % (upi, taxid),
            upi,
            taxid,
            int(self.has_issue),
            int(self.incomplete_sequence),
            int(self.possible_contamination),
            int(self.missing_rfam_match),
            psql.list_as_array(self.messages),
        ]
