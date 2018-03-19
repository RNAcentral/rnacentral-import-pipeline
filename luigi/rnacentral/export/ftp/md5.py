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

from rnacentral.psql import PsqlWrapper

KNOWN = """
select upi, md5 from rna order by id asc
"""


def known(config, handle):
    """
    Generate a file of the UPI to MD5 mappings.
    """

    psql = PsqlWrapper(config)
    psql.write_query(handle, KNOWN, use='tsv')


def example(config, handle, count):
    """
    Generate the example file of UPI to MD5 mappings. This will given count
    number of the first mappings.
    """

    psql = PsqlWrapper(config)
    psql.write_query(handle, KNOWN + " limit %s" % count, use='tsv')
