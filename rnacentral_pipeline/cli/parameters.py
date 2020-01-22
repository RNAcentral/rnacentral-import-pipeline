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

import click

from rnacentral_pipeline.rnacentral.traveler import data as traveler


class TravelerSourceParam(click.ParamType):
    name = 'traveler-source'

    def convert(self, value, param, ctx):
        if value is None:
            return None

        if isinstance(value, traveler.Source):
            return value

        if not isinstance(value, str):
            self.fail("Cannot convert: %s" % value)

        name = value.lower()
        if not hasattr(traveler.Source, name):
            self.fail("Unknown source: " + value)
        return getattr(traveler.Source, name)
