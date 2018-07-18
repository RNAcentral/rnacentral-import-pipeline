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

from rnacentral_pipeline.db import cursor


def unique_sequences(db_url):
    with cursor(db_url) as cur:
        cur.execute("select count(distinct upi) from xref where deleted = 'N'")
        return cur.fetchone()[0]


def active_xrefs(db_url):
    with cursor(db_url) as cur:
        cur.execute("select count(*) from xref where deleted = 'N'")
        return cur.fetchone()[0]


def write(template_file, release, output, db_url):
    template = template_file.read()
    output.write(template.format(
        release=release,
        active_xrefs=active_xrefs(db_url),
        unique_sequences=unique_sequences(db_url),
    ))
