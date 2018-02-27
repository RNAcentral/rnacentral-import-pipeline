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

from tasks.config import ena

from tasks.utils.fetch import FetchTask


def tpa_tasks():
    """
    Create a list of tasks to fetch TPA files.
    """

    tasks = []
    conf = ena()
    for database in ena().tpa_databases:
        local_path = conf.input_tpa_file(database)
        remote_path = conf.raw_tpa_url(database)
        tasks.append(FetchTask(remote_path=remote_path, local_path=local_path))
    return tasks


def copy_ncr_task(filename):
    """
    This will create a task to fetch either an NCR file or a directory of NCR
    files.
    """

    conf = ena()
    head, tail = os.path.split(filename)
    _, first = os.path.split(head)
    local_dir = conf.input_ncr_file(first, tail)
    return FetchTask(remote_path=filename, local_path=local_dir)
