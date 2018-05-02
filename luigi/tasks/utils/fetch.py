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
import shutil

import luigi

from . import http
from . import ftp
from . import mysql


def fetch(remote_path, local_path):
    try:
        os.makedirs(os.path.dirname(local_path))
    except:
        pass

    if remote_path.startswith('http'):
        http.download(remote_path, local_path)

    elif remote_path.startswith('ftp:'):
        ftp.download(remote_path, local_path)

    elif remote_path.startswith('mysql:'):
        mysql.download(remote_path, local_path)

    elif os.path.exists(remote_path):
        if os.path.isdir(remote_path):
            shutil.copytree(remote_path, local_path)
        elif os.path.isfile(remote_path):
            shutil.copy(remote_path, local_path)
        else:
            raise ValueError("Unknown type of file")

    else:
        raise ValueError("Remote must be http/ftp url or existing path")


class FetchTask(luigi.Task):
    """
    A class to handle the generic task of downloading a file. The file can be
    accessed via http, ftp or
    """

    remote_path = luigi.Parameter()
    local_path = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.local_path)

    def run(self):
        return fetch(self.remote_path, self.local_path)
