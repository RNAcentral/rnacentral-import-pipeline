# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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


@attr.s()
class RfamClan(object):
    id = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str))
    description = attr.ib(validator=is_a(str))
    families = attr.ib(validator=is_a(set))

    @property
    def family_count(self):
        return len(self.families)


def load_clans(version='CURRENT'):
    clan_file = get_clans(version=version)
    membership_file = get_clan_membership(version=version)
    return RfamClan.build_all(clan_file, membership_file)
