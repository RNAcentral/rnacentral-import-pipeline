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
from contextlib import contextmanager

import pytest

from click.testing import CliRunner

from rnacentral_pipeline.cli import europepmc


@contextmanager
def run_indexed_command(cmd, **kwargs):
    runner = CliRunner()
    directory = os.path.abspath('data/publications/')
    with runner.isolated_filesystem():
        result = runner.invoke(europepmc.cli, ['index-xml', directory])
        assert result.exit_code == 0
        result = runner.invoke(europepmc.cli, cmd, **kwargs)
        yield result


def test_can_index_directory_of_xml():
    runner = CliRunner()
    directory = os.path.abspath('data/publications/')
    with runner.isolated_filesystem():
        result = runner.invoke(europepmc.cli, ['index-xml', directory])
        assert result.exit_code == 0
        assert os.path.exists('references.db')


def test_can_do_lookup_with_fallback():
    ids = [
        'PMID:1903816,AAAA,BBBB',
        'PMID:375006,FIRST,SECOND',
    ]
    cmd = ['lookup', '--allow-fallback', 'references.db', '-']
    with run_indexed_command(cmd, input='\n'.join(ids)) as result:
        print(result.output)
        assert result.exit_code == 0
        assert os.path.exists('references.csv')
        with open('references.csv') as raw:
            assert raw.readlines() == [
                '42555a855f56f8d5001ee07d95fdcf46,AAAA,BBBB,"Leung J, Sinclair DA, Hayashi S, Tener GM, Grigliatti TA.",J Mol Biol 219(2):175-188 (1991),Informational redundancy of tRNA(4Ser) and tRNA(7Ser) genes in Drosophila melanogaster and evidence for intergenic recombination,1903816,10.1016/0022-2836(91)90560-s\n',
                'da0c9805cab7efd21a0eed0e17a8223a,FIRST,SECOND,"Macino G, Tzagoloff A.",Mol Gen Genet 169(2):183-188 (1979),Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae,375006,10.1007/bf00271669\n',
            ]


def test_can_use_specified_column():
    ids = ['FIRST,PMID:375006,SECOND']
    cmd = ['lookup', '--allow-fallback', '--column=1', 'references.db', '-']
    with run_indexed_command(cmd, input='\n'.join(ids)) as result:
        assert result.exit_code == 0
        assert os.path.exists('references.csv')
        with open('references.csv') as raw:
            assert raw.readlines() == [
                'da0c9805cab7efd21a0eed0e17a8223a,FIRST,SECOND,"Macino G, Tzagoloff A.",Mol Gen Genet 169(2):183-188 (1979),Assembly of the mitochondrial membrane system: two separate genes coding for threonyl-tRNA in the mitochondrial DNA of Saccharomyces cerevisiae,375006,10.1007/bf00271669\n'
            ]


def test_will_fail_with_badly_formmated_ref_id():
    ids = ['b,PMID:1903816']
    cmd = ['lookup', '--allow-fallback', 'references.db', '-']
    with run_indexed_command(cmd, input='\n'.join(ids)) as result:
        assert result.exit_code != 0
        assert os.path.exists('references.csv')


def test_will_ignore_lookup_failures():
    ids = [ 'PMID:26184978,b', 'PMID:1903816,a']
    cmd = ['lookup', '--no-allow-fallback', 'references.db', '-']
    with run_indexed_command(cmd, input='\n'.join(ids)) as result:
        assert result.exit_code == 0
        assert os.path.exists('references.csv')


def test_will_fail_lookup_faiures_if_requested():
    ids = ['PMID:26184978,b', 'PMID:1903816,a']
    cmd = ['lookup', '--no-allow-fallback', '--no-ignore-missing', 'references.db', '-']
    with run_indexed_command(cmd, input='\n'.join(ids)) as result:
        assert result.exit_code == -1
        assert os.path.exists('references.csv')
