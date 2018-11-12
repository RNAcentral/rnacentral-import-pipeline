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

# Fetch Ensembl
bin/fetch generic 'ftp://ftp.ensembl.org/pub/current_embl/homo_sapiens/Homo_sapiens.GRCh38.*.chromosome.12.dat.gz' 'data/ensembl/Homo_sapiens.GRCh38.chromosome.12.dat.gz'
bin/fetch generic 'ftp://ftp.ensembl.org/pub/current_embl/homo_sapiens/Homo_sapiens.GRCh38.*.chromosome.X.dat.gz' 'data/ensembl/Homo_sapiens.GRCh38.chromosome.X.dat.gz'
bin/fetch generic 'ftp://ftp.ensembl.org/pub/current_embl/mus_musculus/Mus_musculus.GRCm38.*.chromosome.3.dat.gz' 'data/ensembl/Mus_musculus.GRCm38.chromosome.3.dat.gz'
bin/fetch generic 'ftp://ftp.ensembl.org/pub/current_embl/macaca_mulatta/Macaca_mulatta.Mmul_8.0.1.*.chromosome.1.dat.gz' 'data/ensembl/Macaca_mulatta.Mmul_8.0.1.chromosome.1.dat.gz'
gzip -fd data/ensembl/*.gz

# Fetch GENCODE
gencode_urls="$(mktemp)"
cat >$gencode_urls <<EOF
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz
EOF

bin/fetch gencode "$gencode_urls" 'data/gencode/human-transcripts.gff3'

# Fetch EnsemblPlants
bin/fetch generic 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.*.chromosome.2.dat.gz' 'data/ensembl_plants/Arabidopsis_thaliana.TAIR10.chromosome.2.dat.gz'
bin/fetch generic 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/hordeum_vulgare/Hordeum_vulgare.IBSC_v2.*.chromosome.Pt.dat.gz' 'data/ensembl_plants/Hordeum_vulgare.IBSC_v2.chromosome.Pt.dat.gz'
bin/fetch generic 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/oryza_barthii/Oryza_barthii.O.barthii_v1.*.chromosome.9.dat.gz' 'data/ensembl_plants/Oryza_barthii.O.barthii_v1.chromosome.9.dat.gz'
bin/fetch generic 'ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/zea_mays/Zea_mays.B73_RefGen_v4.*.chromosome.7.dat.gz' 'data/ensembl_plants/Zea_mays.B73_RefGen_v4.chromosome.7.dat.gz'
gzip -fd data/ensembl_plants/*.gz
