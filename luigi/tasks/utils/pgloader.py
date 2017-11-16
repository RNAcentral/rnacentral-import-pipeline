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
import abc
from hashlib import sha256

from luigi.local_target import LocalTarget
from luigi.local_target import atomic_file
from luigi.contrib.external_program import ExternalProgramTask

from utils import snake_case
from tasks.config import output
from tasks.config import db as DBConfig


class PGLoader(ExternalProgramTask):  # pylint: disable=R0921
    """
    This is a base class for other classes to run pgloader. The subclass must
    implement the control_file method which generates the text of the pgloader
    control file. In addition, the class must have a unique name, as the
    control file is named using it. The presence or absence of this control
    file is what determines if the task is complete or not.

    PGloader will log to a file called logs/class_name.log, relative to the
    globally configured base directory (see config.output). This log file will
    be appended to each run, not overwritten. It is also ignored when
    determining if the task is complete.
    """
    __metaclass__ = abc.ABCMeta

    def control_filename(self):
        """
        Create the filename for the control file. The file will be named using
        a snake_case version of the class name, and placed in the cmds
        subdirectory of the configured base directory.
        """
        return self.__directory_filename__('cmds')

    def output(self):
        filename = self.control_filename()
        dirname = os.path.dirname(filename)
        try:
            os.makedirs(dirname)
        except Exception as err:
            if not os.path.exists(dirname):
                raise err
        return LocalTarget(filename)

    def program_args(self):
        log_file = self.__directory_filename__('logs', suffix='log')
        return ['pgloader', '--logfile', log_file, self.control_filename()]

    @abc.abstractmethod
    def control_file(self):
        """
        Generate the text of the control file. The `db_url` will be filled out
        prior to writing the text.

        :returns str: The text to format.
        """
        return ''

    def db_url(self, table=None):
        """
        Create a pgloader url for the postgres database that is configured in
        the db section. If the table is given and not None then it will be
        added as well.
        """

        db_url = DBConfig().pgloader_url()
        if table is not None:
            return '%s?%s' % (db_url, table)
        return db_url

    def db_search_path(self):
        """
        Get the configured search path.
        """
        return DBConfig().search_path

    def write_control_file(self):
        """
        This will write the control file just prior to running pgloader. This
        will also format the text returned by `control_file` to fill out the
        `db_url` parameter with the url for the postgres database as configured
        by config.db.
        """
        filename = self.control_filename()
        with atomic_file(filename) as out:
            out.write(self.control_file())
        return filename

    def run(self):
        """
        This will run pgloader on the generated control file. If pgloader runs
        with no issues, then this will keep the file, otherwise the file will
        be renamed by adding a '-crashed' prefix to assist with debugging.
        """
        filename = self.write_control_file()
        try:
            super(PGLoader, self).run()
        except:
            os.rename(filename, '%s-crashed' % filename)
            raise

    def content_hash(self):
        """
        This produces a hash of the content of the control file. This is useful
        in naming the files as sometimes the file depends on some given
        parameters and instead of forcing the user to override the filename
        written so we always have a unique filename this will provide a unique
        hash for each parameter.
        """

        content = self.control_file()
        return sha256(content).hexdigest()

    def __directory_filename__(self, directory, suffix='ctl'):
        directory = os.path.join(output().base, directory)
        name = snake_case(self.__class__.__name__)
        filename = '%{name}-{hash}.{suffix}'.format(
            name=name,
            hash=self.content_hash(),
            suffix=suffix,
        )
        return os.path.join(directory, filename)
