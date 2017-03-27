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

from parameters import GenericFileParameter


class Download(luigi.Task):
    """
    A simple class to download remote files and store them locally. They will
    """
    remote_file = GenericFileParameter()
    directory = luigi.Parameter(default='./data')

    @property
    def remote(self):
        return GenericFileParameter().as_target(self.remote_file)

    def output(self):
        return LocalTarget(path=self.path())

    def path(self):
        remote = self.remote.path
        if remote.endswith('.gz'):
            remote = remote[:-3]
        if remote.startswith('/'):
            remote = remote[1:]
        return os.path.join(self.directory, remote)

    def run(self):
        path = self.path()
        basedir = os.path.dirname(path)
        if not os.path.exists(basedir):
            os.makedirs(basedir)

        with self.remote.open('r') as raw, atomic_file(path) as out:
            out.writelines(r for r in raw)


if __name__ == '__main__':
    luigi.run(main_task_cls=Download)
