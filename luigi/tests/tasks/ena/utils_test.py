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

from tasks.config import output
from tasks.ena import utils


def test_tpa_task_are_generated_correctly():
    tasks = utils.tpa_tasks()
    assert len(tasks) == 10
    tasks[0].remote_path == 'https://www.ebi.ac.uk/ena/data/xref/search?source=PomBase&expanded=true'
    tasks[0].local_path == os.path.join(output().base, 'ena', 'tpa', 'PomBase.tsv')


def test_can_produce_copy_for_nested_directories_tasks():
    task = utils.copy_ncr_task('/nfs/ftp/pub/databases/ena/non-coding/release/wgs/aj')
    assert task.local_path == os.path.join(output().base, 'ena', 'ncr', 'wgs', 'aj')
    assert task.remote_path == '/nfs/ftp/pub/databases/ena/non-coding/release/wgs/aj'


def test_can_produce_copy_for_simple_directories():
    task = utils.copy_ncr_task('/nfs/ftp/pub/databases/ena/non-coding/release/std/rel_est_rod_03_r134.ncr.gz')
    assert task.local_path == os.path.join(output().base, 'ena', 'ncr', 'std', 'rel_est_rod_03_r134.ncr.gz')
    assert task.remote_path == '/nfs/ftp/pub/databases/ena/non-coding/release/std/rel_est_rod_03_r134.ncr.gz'
