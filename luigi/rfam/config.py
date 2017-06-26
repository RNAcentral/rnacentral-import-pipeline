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

import luigi
from parameters import PathParameter


class output(luigi.Config):
    base = PathParameter(default='/tmp')


class db(luigi.Config):
    """
    This contains the configuration information about the database to connect
    to.
    """

    user = luigi.Parameter(default='rnacen')
    password = luigi.Parameter()
    host = luigi.Parameter(default='127.0.0.1')
    port = luigi.Parameter(default=5432)
    db_name = luigi.Parameter(default='rnacentral')


class files(luigi.Config):
    """
    This contains the configuration about the rfam files to read. Notably, the
    rfam_hits value is the path to the .tbl file from an Rfam search to import.
    """
    rfam_hits = PathParameter()
