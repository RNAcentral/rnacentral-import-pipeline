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

import csv
import itertools as it

from . import databases


def load_done(handle):
    return {tuple(r) for r in csv.reader(handle)}


def load_possible(handle):
    return {line.strip() for line in handle}


def find_tasks(handle, possible_handle, done_handle):
    done = load_done(done_handle)
    possible = load_possible(possible_handle)
    dbs = list(databases.select_max(handle))
    todo = it.product(possible, dbs)
    return [t for t in todo if t not in done]


def write(handle, possible_handle, done_handle, output):
    tasks = find_tasks(handle, possible_handle, done_handle)
    csv.writer(output).writerows(tasks)
