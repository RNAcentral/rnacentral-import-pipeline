# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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


@attr.s()
class Interaction:
    external_id: str = attr.ib(validator=is_a(str))
    source: str = attr.ib(validator=is_a(str))


    def writeables(self):
        yield [
            self.external_id,
            self.source,
        ]
