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

import json

import semver

from . import v1

from rnacentral_pipeline.writers import build_entry_writer


def parse(handle):
    """
    This parses the file like object that should contain the RNAcentral data.
    The file can contain data in any of the accepted versions.
    """

    data = json.load(handle)
    if not data:
        raise ValueError("No data loaded")
    if not data['data']:
        raise ValueError("Missing data to import")

    version = data.get('metaData', {}).get('schemaVersion', None)
    if not version:
        raise ValueError("Must specify a schema version in metadata")

    if semver.match(version, "<=2.0.0"):
        return v1.parse(data)
    raise ValueError("Unknown schema version: %s" % version)


def from_file(handle, output):
    writer = build_entry_writer(parse)
    writer(output, handle)
