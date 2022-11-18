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

import os
import random
from contextlib import contextmanager

import pytest
from click.testing import CliRunner
from more_itertools import chunked

from rnacentral_pipeline.cli import genome_mapping as gm


@contextmanager
def empty_input(args, **kwargs):
    runner = CliRunner()
    with runner.isolated_filesystem():
        open("empty", "w").close()
        assert os.path.exists("empty")
        with open("empty", "r") as raw:
            assert raw.read() == ""
        result = runner.invoke(gm.cli, args, **kwargs)
        yield result


@pytest.mark.cli
@pytest.mark.parametrize(
    "command,extra",
    [
        ("serialize", "a"),
        ("as-importable", None),
        ("select", None),
    ],
)
def test_blat_commands_handle_empty_inputs(command, extra):
    cmd = ["blat", command, "empty", "output"]
    if extra:
        cmd.insert(2, extra)
    with empty_input(cmd) as result:
        assert result.exit_code == 0
        assert os.path.exists("output")
        with open("output", "r") as raw:
            assert raw.read() == ""


@pytest.mark.cli
@pytest.mark.parametrize(
    "filename",
    [
        os.path.join("data/genome-mapping/", f)
        for f in os.listdir("data/genome-mapping/")
    ],
)
def test_always_produces_output(filename):
    runner = CliRunner()
    filename = os.path.abspath(filename)
    with runner.isolated_filesystem():
        assert os.path.exists(filename)
        cmd = ["blat", "serialize", "a", filename, "output.pickle"]
        result = runner.invoke(gm.cli, cmd)
        assert result.exit_code == 0
        assert os.path.exists("output.pickle")

        cmd = ["blat", "select", "output.pickle", "selected.pickle"]
        runner.invoke(gm.cli, cmd)
        assert result.exit_code == 0
        assert os.path.exists("selected.pickle")

        cmd = ["blat", "as-importable", "selected.pickle", "selected.csv"]
        runner.invoke(gm.cli, cmd)
        assert result.exit_code == 0
        assert os.path.exists("selected.csv")


@pytest.mark.cli
def test_sorting_works_correctly():
    runner = CliRunner()
    filename = os.path.abspath("data/genome-mapping/results.psl")
    parts = 1000
    with open(filename, "r") as raw:
        lines = raw.readline()
        lines = random.shuffle(lines)

    with runner.isolated_filesystem():
        for index, chunk in enumerate(chunked(lines, parts)):
            chunk_name = "results-%i" % index
            psl = "%s.psl"
            stored = "%s.pickle"
            selected = "%s-selected.pickle"
            with open(psl, "w") as out:
                out.writelines(chunk)

            cmd = ["blat", "serialize", chunk_name, stored]
            runner.invoke(gm.cli, cmd)
            assert result.exit_code == 0
            assert os.path.exists(stored)
