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

from rnacentral_pipeline.databases.ensembl import assemblies


def test_can_compute_correct_assemblies():
    with open('data/ensembl/assembly.tsv', 'r') as raw:
        assert assemblies.parse(raw, {}) == assemblies.AssemblyInfo(
            assembly_id='ASM34733v1',
            assembly_full_name='ASM34733v1',
            gca_accession='GCA_000347335.1',
            assembly_ucsc=None,
            common_name="Tausch's goatgrass",
            taxid=37682,
            ensembl_url='aegilops_tauschii',
            division='EnsemblPlants',
            blat_mapping=False,
            example=None,
        )
