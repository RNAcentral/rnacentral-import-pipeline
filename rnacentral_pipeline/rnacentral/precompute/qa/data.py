# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import json
import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


@attr.s(frozen=True)
class QaResult:
    """
    This represents the result of a single validator.
    """

    name = attr.ib(validator=is_a(str))
    has_issue = attr.ib(validator=is_a(bool))
    message = attr.ib(validator=optional(is_a(str)))

    @classmethod
    def ok(cls, name):
        return cls(name=name, has_issue=False, message=None)

    @classmethod
    def not_ok(cls, name, message):
        return cls(name=name, has_issue=True, message=message)

    def str_issue(self):
        return str(int(self.has_issue))


@attr.s(frozen=True)
class QaStatus:
    """
    This represents an update to the QA status table.
    """

    incomplete_sequence = attr.ib(validator=is_a(QaResult))
    possible_contamination = attr.ib(validator=is_a(QaResult))
    missing_rfam_match = attr.ib(validator=is_a(QaResult))
    from_repetitive_region = attr.ib(validator=is_a(QaResult))

    @classmethod
    def from_results(cls, results: ty.List[QaResult]) -> "QaStatus":
        fields = attr.fields_dict(cls)
        data = {}
        for result in results:
            if result.name not in fields:
                raise ValueError(f"Unknown QaResult {result}")
            data[result.name] = result
        if not data:
            raise ValueError("Cannot build without QA results")
        return cls(**data)

    @property
    def has_issue(self) -> bool:
        """
        Check if this QA update indicates if there is any issue.
        """

        return (
            self.incomplete_sequence.has_issue
            or self.possible_contamination.has_issue
            or self.missing_rfam_match.has_issue
            or self.from_repetitive_region.has_issue
        )

    def messages(self) -> ty.List[str]:
        messages = []
        for field in attr.fields(self.__class__):
            result = getattr(self, field.name)
            if result.message:
                messages.append(result.message)
        return messages

    def writeable(self, upi: str, taxid: int) -> ty.List[str]:
        """
        Create a writeable array for writing CSV files.
        """

        return [
            "%s_%i" % (upi, taxid),
            upi,
            str(taxid),
            str(int(self.has_issue)),
            self.incomplete_sequence.str_issue(),
            self.possible_contamination.str_issue(),
            self.missing_rfam_match.str_issue(),
            json.dumps(self.messages()),
        ]
