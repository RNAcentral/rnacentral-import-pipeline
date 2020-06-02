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

import re
import itertools as it

from Bio import SeqIO

NHMMER_PATTERN = re.compile('^[ABCDGHKMNRSTVWXYU]+$', re.IGNORECASE)


def is_valid_nhmmer_record(record):
    """
    Checks if a sequence is valid for nhmmer usage.
    """
    return bool(NHMMER_PATTERN.match(str(record.seq)))


def valid_nhmmer(handle, output):
    sequences = SeqIO.parse(handle, 'fasta')
    accepted = filter(is_valid_nhmmer_record, sequences)
    SeqIO.write(accepted, output, 'fasta')


def invalid_nhmmer(handle, output):
    sequences = SeqIO.parse(handle, 'fasta')
    rejected = it.filterfalse(is_valid_nhmmer_record, sequences)
    SeqIO.write(rejected, output, 'fasta')
