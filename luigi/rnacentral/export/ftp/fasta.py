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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rnacentral.psql import PsqlWrapper

BASE_ACTIVE = """
select
    pre.id,
    pre.description,
    case when rna.seq_short is null
        then rna.seq_long
        else rna.seq_short
    end as sequence
from rna_active active
join rna on rna.upi = active.upi
join rnc_rna_precomputed pre
on pre.upi = rna.upi and pre.upi = active.upi
"""

ACTIVE_SQL = BASE_ACTIVE + " where pre.taxid is null"


EXAMPLE_SQL = BASE_ACTIVE + " order by pre.upi limit 10"


ACITVE_SPECIES_SQL = BASE_ACTIVE + " where pre.taxid is not null"


INACTIVE_SQL = """
select
    pre.id,
    pre.description,
    case when rna.seq_short is null
        then rna.seq_long
        else rna.seq_short end as sequence
from rna_inactive inactive
join rna on rna.upi = inactive.upi
join rnc_rna_precomputed pre
on pre.upi = rna.upi and pre.upi = inactive.upi
where
    pre.taxid is null
"""


NHMMER_PATTERN = re.compile('^[ABCDGHKMNRSTVWXYU]+$', re.IGNORECASE)


def as_record(entry):
    return SeqRecord(
        Seq(entry['sequence']),
        id=entry['id'],
        description=entry['description'],
    )


def is_valid_nhmmer_record(record):
    return bool(NHMMER_PATTERN.match(str(record.seq)))


def active(config):
    psql = PsqlWrapper(config)
    sql = ACTIVE_SQL + " order by pre.upi"
    return it.imap(as_record, psql.copy_to_iterable(sql))


def species(config):
    psql = PsqlWrapper(config)
    sql = ACITVE_SPECIES_SQL + " order by pre.upi"
    return it.imap(as_record, psql.copy_to_iterable(sql))


def example(config):
    psql = PsqlWrapper(config)
    return it.imap(as_record, psql.copy_to_iterable(EXAMPLE_SQL))


def inactive(config):
    psql = PsqlWrapper(config)
    sql = INACTIVE_SQL + " order by pre.upi"
    return it.imap(as_record, psql.copy_to_iterable(sql))


def nhmmer_records(filename, select_valid=True):
    for record in SeqIO.parse(filename, 'fasta'):
        if select_valid == is_valid_nhmmer_record(record):
            yield record
