# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import pytest

from rnacentral_pipeline.databases.ensembl.vertebrates import urls


class FakeFtp:
    """
    Only serves the paths it is given, so asking for anything else raises the
    way the real server does with a 550.
    """

    def __init__(self, files):
        self.files = files

    def retrlines(self, cmd, callback):
        _, path = cmd.split(" ", 1)
        if path not in self.files:
            raise Exception(f"550 Failed to open file: {path}")
        for line in self.files[path].splitlines():
            callback(line)


def test_reads_release_from_version_not_readme():
    """
    current_README was removed from ftp.ensembl.org; pub/VERSION holds the
    release number now.
    """
    assert urls.latest_release(FakeFtp({"VERSION": "116\n"})) == "release-116"


def test_rejects_a_version_that_is_not_a_release_number():
    ftp = FakeFtp({"VERSION": "<html>404 Not Found</html>"})
    with pytest.raises(ValueError):
        urls.latest_release(ftp)
