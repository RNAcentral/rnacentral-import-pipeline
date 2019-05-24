#!/usr/bin/env bash

# Copyright [2009-2017] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -euo pipefail
IFS=$'\n\t'

fetch()
{
  echo $2
  wget -q -O - $1 | gzip -d - > $2
}

fetch 'ftp://ftp.ensembl.org/pub/current_embl/homo_sapiens/Homo_sapiens.GRCh38.*.chromosome.3.dat.gz' 'data/ensembl/Homo_sapiens.GRCh38.chromosome.3.dat'
fetch 'ftp://ftp.ensembl.org/pub/current_embl/homo_sapiens/Homo_sapiens.GRCh38.*.chromosome.12.dat.gz' 'data/ensembl/Homo_sapiens.GRCh38.chromosome.12.dat'
fetch 'ftp://ftp.ensembl.org/pub/current_embl/homo_sapiens/Homo_sapiens.GRCh38.*.chromosome.X.dat.gz' 'data/ensembl/Homo_sapiens.GRCh38.chromosome.X.dat'
fetch 'ftp://ftp.ensembl.org/pub/current_embl/mus_musculus/Mus_musculus.GRCm38.*.chromosome.3.dat.gz' 'data/ensembl/Mus_musculus.GRCm38.chromosome.3.dat'
fetch 'ftp://ftp.ensembl.org/pub/current_embl/macaca_mulatta/Macaca_mulatta.Mmul_8.0.1.*.chromosome.1.dat.gz' 'data/ensembl/Macaca_mulatta.Mmul_8.0.1.chromosome.1.dat'
fetch 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz' 'data/gencode/human-transcripts.gff3'
fetch 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.*.chromosome.2.dat.gz' 'data/ensembl_plants/Arabidopsis_thaliana.TAIR10.chromosome.2.dat'
fetch 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/hordeum_vulgare/Hordeum_vulgare.IBSC_v2.*.chromosome.Pt.dat.gz' 'data/ensembl_plants/Hordeum_vulgare.IBSC_v2.chromosome.Pt.dat'
fetch 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/oryza_barthii/Oryza_barthii.O.barthii_v1.*.chromosome.9.dat.gz' 'data/ensembl_plants/Oryza_barthii.O.barthii_v1.chromosome.9.dat'
fetch 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/zea_mays/Zea_mays.B73_RefGen_v4.*.chromosome.7.dat.gz' 'data/ensembl_plants/Zea_mays.B73_RefGen_v4.chromosome.7.dat'
fetch 'ftp://ftp.ensembl.org/pub/release-96/embl/bos_taurus/Bos_taurus.ARS-UCD1.2.96.primary_assembly.8.dat.gz' 'data/ensembl/Bos_taurus.ARS-UCD1.2.primary_assembly.8.dat'
