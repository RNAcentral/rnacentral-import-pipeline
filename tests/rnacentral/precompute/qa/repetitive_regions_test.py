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

from types import SimpleNamespace

from rnacentral_pipeline.rnacentral.precompute.qa import repetitive_regions


def coordinate(assembly_id="GRCh38", chromosome="1", start=100, stop=200):
    return SimpleNamespace(
        assembly_id=assembly_id, chromosome=chromosome, start=start, stop=stop
    )


class FakeAssemblyRepeats:
    def __init__(self, enveloped):
        self.enveloped = enveloped

    def is_enveloped(self, chromosome, start, stop):
        return self.enveloped


class FakeRepeats:
    def __init__(self, assemblies=None):
        self.assemblies = assemblies or {}

    def has_assembly(self, assembly_id):
        return assembly_id in self.assemblies

    def assembly(self, assembly_id):
        return self.assemblies[assembly_id]


def context(assemblies=None):
    return SimpleNamespace(repeats=FakeRepeats(assemblies))


def sequence(*coordinates):
    return SimpleNamespace(coordinates=list(coordinates))


def test_ok_when_sequence_has_no_coordinates():
    result = repetitive_regions.validate(context(), sequence())
    assert result.has_issue is False


def test_ok_when_assembly_has_no_repeat_data():
    ctx = context({})
    result = repetitive_regions.validate(ctx, sequence(coordinate()))
    assert result.has_issue is False


def test_ok_when_coordinate_is_not_enveloped_by_a_repeat():
    ctx = context({"GRCh38": FakeAssemblyRepeats(False)})
    result = repetitive_regions.validate(ctx, sequence(coordinate()))
    assert result.has_issue is False


def test_flags_issue_when_coordinate_is_enveloped_by_a_repeat():
    ctx = context({"GRCh38": FakeAssemblyRepeats(True)})
    result = repetitive_regions.validate(ctx, sequence(coordinate()))
    assert result.has_issue is True
    assert result.message == "This sequence overlaps a repetitive region"


def test_flags_issue_if_any_coordinate_is_enveloped():
    ctx = context(
        {
            "GRCh38": FakeAssemblyRepeats(False),
            "GRCh37": FakeAssemblyRepeats(True),
        }
    )
    coords = sequence(
        coordinate(assembly_id="GRCh38"),
        coordinate(assembly_id="GRCh37"),
    )
    result = repetitive_regions.validate(ctx, coords)
    assert result.has_issue is True
