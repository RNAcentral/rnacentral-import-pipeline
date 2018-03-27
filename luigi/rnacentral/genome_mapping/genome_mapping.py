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

import csv
import subprocess
import sys

from rnacentral.psql import PsqlWrapper

import itertools as it

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



FASTA = """
with ensembl_ids as (
    select distinct upi
    from xref
    where dbid = 25
    and deleted = 'N'
    and taxid = {taxid}
), no_ensembl as (
    select distinct xref.upi, xref.taxid
    from xref
    left outer join ensembl_ids
    on xref.upi = ensembl_ids.upi
    where xref.deleted = 'N'
    and xref.taxid = {taxid}
    and ensembl_ids.upi is null
)
select
	concat_ws('_', t1.upi, t2.taxid) as id,
	(CASE WHEN (t1.seq_long IS NULL) THEN t1.seq_short ELSE t1.seq_long END) as sequence
from rna t1, no_ensembl t2
where
t1.upi = t2.upi
"""

def as_record(entry):
    """
    Turns the entry into a SeqRecord for output.
    """
    return SeqRecord(
        Seq(entry['sequence']).transcribe(),
        id=entry['id'],
        description='',
    )

def export_rnacentral_fasta(config, taxid):
    """
    Return an iteratable of sequence records.
    """
    psql = PsqlWrapper(config)
    sql = FASTA.format(taxid=taxid)
    return it.imap(as_record, psql.copy_to_iterable(sql))
