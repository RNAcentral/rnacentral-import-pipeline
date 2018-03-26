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
from luigi import LocalTarget
from luigi.local_target import atomic_file

from tasks.config import output
from tasks.utils.parameters import GenericFileParameter


class Download(luigi.Task):  # pylint: disable=R0904
    """
    A simple class to download remote files and store them locally. They will
    """

    remote_file = GenericFileParameter()

    @property
    def remote(self):
        """
        Create a Target for the remote file.
        """
        return GenericFileParameter().as_target(self.remote_file)

    def output(self):
        return LocalTarget(path=self.path())

    def path(self):
        """
        Determine the local path to write to.
        """

        remote = self.remote.path
        if remote.endswith('.gz'):
            remote = remote[:-3]
        if remote.startswith('/'):
            remote = remote[1:]
        return os.path.abspath(os.path.join(output().base, remote))

    def run(self):
        """
        Download the file. This will write the file, using the same path
        structure as the remote location has.
        """

        path = self.path()
        basedir = os.path.dirname(path)
        try:
            os.makedirs(basedir)
        except:
            pass

        with self.remote.open('r') as raw, atomic_file(path) as out:
            out.writelines(r for r in raw)
