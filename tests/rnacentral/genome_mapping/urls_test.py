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

import six

import pytest

from rnacentral_pipeline.rnacentral.genome_mapping import urls


@pytest.mark.parametrize('name,expected', [
    ('ensembl', urls.FtpHost.ensembl),
    ('ensembl_vertebrates', urls.FtpHost.ensembl),
    ('ensembl_genomes', urls.FtpHost.ensembl_genomes),
    ('ensembl_plants', urls.FtpHost.ensembl_plants),
    ('EnsemblPlants', urls.FtpHost.ensembl_plants),
    ('unknown', urls.FtpHost.unknown),
])
def test_can_generate_hosts_correctly(name, expected):
    assert urls.FtpHost.from_string(name) is expected


@pytest.mark.parametrize('species,assembly_id,host,expected', [
    ('amphiprion_ocellaris', 'AmpOce1.0', urls.FtpHost.ensembl,
     'ftp://ftp.ensembl.org/pub/current_fasta/amphiprion_ocellaris/dna/Amphiprion_ocellaris.AmpOce1.0.dna.toplevel.fa.gz'),
    ('anolis_carolinensis', 'AnoCar2.0', urls.FtpHost.ensembl,
     'ftp://ftp.ensembl.org/pub/current_fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa.gz'),
    ('homo_sapiens', 'GRCh38', urls.FtpHost.ensembl,
     'ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'),
    ('drosophila_sechellia', 'dsec_caf1', urls.FtpHost.ensembl_metazoa,
     'ftp://ftp.ensemblgenomes.org/pub/current/metazoa/fasta/drosophila_sechellia/dna/Drosophila_sechellia.dsec_caf1.dna.toplevel.fa.gz'),
    ('drosophila_sechellia', 'dsec_caf1', urls.FtpHost.ensembl_genomes,
     'ftp://ftp.ensemblgenomes.org/pub/current/metazoa/fasta/drosophila_sechellia/dna/Drosophila_sechellia.dsec_caf1.dna.toplevel.fa.gz'),
    ('drosophila_sechellia', 'dsec_caf1', urls.FtpHost.unknown,
     'ftp://ftp.ensemblgenomes.org/pub/current/metazoa/fasta/drosophila_sechellia/dna/Drosophila_sechellia.dsec_caf1.dna.toplevel.fa.gz'),
])
def test_can_find_correct_url_for_species(species, assembly_id, host, expected):
    assert urls.url_for(species, assembly_id, host=host) == expected


@pytest.mark.parametrize('species,assembly_id', [
    ('anas_platyrhynchos', 'CAU_duck1.0'),  # Different species on FTP only
])
def test_raises_exceptions_for_weird_url_cases(species, assembly_id):
    with pytest.raises(urls.NoTopLevelFiles):
        urls.url_for(species, assembly_id)


def test_can_get_all_paths_for_a_csv():
    data = six.StringIO()
    data.write('drosophila_sechellia,dsec_caf1,unknown\n')
    data.write('homo_sapiens,GRCh38,ensembl\n')
    data.write('oryza_sativa,IRGSP-1.0,EnsemblPlants\n')
    data.seek(0)
    assert list(urls.urls_for(data)) == [
        ('drosophila_sechellia', 'dsec_caf1', 'ftp://ftp.ensemblgenomes.org/pub/current/metazoa/fasta/drosophila_sechellia/dna/Drosophila_sechellia.dsec_caf1.dna.toplevel.fa.gz'),
        ('homo_sapiens', 'GRCh38', 'ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'),
        ('oryza_sativa', 'IRGSP-1.0', 'ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz'),
    ]
