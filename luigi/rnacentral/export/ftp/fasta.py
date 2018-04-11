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
with rna_status as (
select
    upi,
    '{{N}}' && array_agg(deleted) is_active
from xref
group by upi
)
select
    pre.id,
    pre.description,
    case when rna.seq_short is null
        then rna.seq_long
        else rna.seq_short
    end as sequence
from rnc_rna_precomputed pre
join rna on rna.upi = pre.upi
join rna_status status on status.upi = rna.upi
where
    status.is_active = {is_active}
    and {terms}
"""

ACTIVE_SQL = BASE_ACTIVE.format(
    terms="pre.taxid is null",
    is_active='true',
)

ACITVE_SPECIES_SQL = BASE_ACTIVE.format(
    terms="pre.taxid is not null",
    is_active='true',
)

INACTIVE_SQL = BASE_ACTIVE.format(
    terms="pre.taxid is null",
    is_active='false',
)

NHMMER_PATTERN = re.compile('^[ABCDGHKMNRSTVWXYU]+$', re.IGNORECASE)


def as_record(entry):
    """
    Turns the entry into a SeqRecord for output.
    """
    return SeqRecord(
        Seq(entry['sequence']).transcribe(),
        id=entry['id'],
        description=entry['description'],
    )


def is_valid_nhmmer_record(record):
    """
    Checks if a sequence is valid for nhmmer usage.
    """
    return bool(NHMMER_PATTERN.match(str(record.seq)))


def active(config):
    """
    Extracts all active sequences and produces an iterable of SeqRecords.
    """

    psql = PsqlWrapper(config)
    sql = ACTIVE_SQL + " order by pre.upi"
    return it.imap(as_record, psql.copy_to_iterable(sql))


def species(config):
    """
    Extracts all active sequences and produces an iterable of SeqRecords with a
    species specific id.
    """

    psql = PsqlWrapper(config)
    sql = ACITVE_SPECIES_SQL + " order by pre.upi"
    return it.imap(as_record, psql.copy_to_iterable(sql))


def inactive(config):
    """
    Extract all inactive sequences and produce an iterable of SeqRecords.
    """

    psql = PsqlWrapper(config)
    sql = INACTIVE_SQL + " order by pre.upi"
    return it.imap(as_record, psql.copy_to_iterable(sql))


def nhmmer(handle, select_valid=True):
    """
    Extract all sequences which may be written to nhmmer from the given file.
    If select_valid is True then only valid sequences are written. If it is
    False then only invalid sequences are written.
    """

    for record in SeqIO.parse(handle, 'fasta'):
        if select_valid == is_valid_nhmmer_record(record):
            yield record
