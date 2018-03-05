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

import requests

from luigi.local_target import atomic_file


def download(url, filename):
    """
    This will fetch some file over HTTP using requests. It will create the
    required directory to save in if requried as well. Note there is a race
    condition in the directory creation, so if that is a problem create it
    ahead of time.
    """

    dirname = os.path.dirname(filename)
    try:
        os.makedirs(dirname)
    except:
        pass

    response = requests.get(url)
    response.raise_for_status()
    with atomic_file(filename) as out:
        for chunk in response.iter_content(chunk_size=128):
            out.write(chunk)
