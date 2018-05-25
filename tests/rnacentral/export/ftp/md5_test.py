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

import tempfile

from rnacentral.export.ftp import md5

from tasks.config import db


def test_can_produce_md5_mappings():
    with tempfile.NamedTemporaryFile() as tmp:
        md5.example(db(), tmp, 5)
        tmp.flush()

        with open(tmp.name, 'r') as raw:
            assert raw.readlines() == [
                'URS0000000001\t6bba097c8c39ed9a0fdf02273ee1c79a\n',
                'URS0000000002\t1fe2f0e3c3a2d6d708ac98e9bfb1d7a8\n',
                'URS0000000003\t7bb11d0572ff6bb42427ce74450ba564\n',
                'URS0000000004\t030c78be0f492872b95219d172e0c658\n',
                'URS0000000005\t030c795b3b5bb84256b0fea3c10ee3aa\n',
            ]
