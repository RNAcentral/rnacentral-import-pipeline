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

import luigi
from luigi.local_target import atomic_file

from tasks.config import export


class FastaReadme(luigi.Task):
    def output(self):
        return luigi.LocalTarget(export().sequences('readme.txt'))

    def run(self):
        try:
            os.makedirs(os.path.dirname(self.output().fn))
        except:
            pass

        with atomic_file(self.output().fn) as out:
            out.write("""
===================================================================
RNAcentral Sequence Data
===================================================================

This directory contains sequences with RNAcentral ids in FASTA format.

* rnacentral_active.fasta.gz
Current set of sequences that are present in at least one expert database.

* rnacentral_species_specific_ids.fasta.gz
Current set of sequences that are present in at least one expert database using
the species specific URS ID's.

* rnacentral_inactive.fasta.gz
All RNAcentral sequences that used to be present in one or more expert database
but don't have any current cross-references.

* example.txt
A small example file showing the format of rnacentral_active.fasta.gz
and rnacentral_inactive.fasta.gz.
            """)
