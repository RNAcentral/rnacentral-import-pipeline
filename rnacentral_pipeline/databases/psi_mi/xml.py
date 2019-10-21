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

from lxml import etree as ET


def parse_2_5(root):
    for interaction in root.findall('interaction'):
        yield interaction


def get_version(root):
    return '2.5'


def parse(handle):
    tree = ET.parse(handle)
    root = tree.getroot()
    version = get_version(root)
    if version == '2.5':
        parse_2_5(root)
    else:
        raise ValueError("Unknown version (%s) of PSI-MI XML" % version)
