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

import os

import luigi

from tasks.config import output
from tasks.utils.pgload_sequences import PGLoadSequences


class LoadShort(PGLoadSequences):  # pylint: disable=R0904
    """
    This will load all short sequences. The database paraemter defaults to all
    sequences, if a value is given then it is assumed to be the name of the
    database to load. All files that begin with that name will be loaded.
    """

    database = luigi.Parameter(default='all')

    def control_filename(self):
        suffix = '%s.ctl' % self.database
        return self.__directory_filename__('cmds', suffix=suffix)

    def pattern(self):
        if self.database == 'all':
            return '.*'
        return self.database + '.*.csv'

    def directory(self):
        return os.path.join(output().base, 'short')
