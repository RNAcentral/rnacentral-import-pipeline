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

from .utils import FastaExportBase


CREATE_ACTIVE_TABLE = """
CREATE TEMP TABLE rna_active (
    upi varchar(30) primary KEY
)
"""

INSERT_INTO_ACTIVE_TABLE = """
INSERT INTO rna_active
select upi from xref_p{dbid}_not_deleted
ON CONFLICT DO NOTHING
"""


class ActiveFastaExport(FastaExportBase):
    """
    This task will create the fasta file of all active RNAcentral sequences. It
    does not create the species specific ids and produces a fasta file with the
    generic UPI's.
    """

    filename = 'rnacentral_active.fasta'

    @property
    def fetch(self):
        return """
        select
            pre.id,
            pre.description,
            case when rna.seq_short is null
                then rna.seq_long
                else rna.seq_short end as sequence
        from rna_active active
        join rnc_rna_precomputed pre on pre.upi = active.upi
        join rna on rna.upi = active.upi
        where
            pre.id = rna.upi
            and pre.upi = active.upi
        """


class SpeciesSpecificFastaExport(FastaExportBase):
    """
    Export all species specific sequences. This produces the
    rnacentral_species_specific_ids.fasta file which contains all active
    sequences with the sequence specific UPI's.
    """

    filename = 'rnacentral_species_specific_ids.fasta'

    @property
    def fetch(self):
        return """
        select
            pre.id,
            pre.description,
            case when rna.seq_short is null
                then rna.seq_long
                else rna.seq_short end as sequence
        from rna_active active
        join rnc_rna_precomputed pre on pre.upi = active.upi
        join rna on active.upi = rna.upi
        where
            pre.taxid is not null
        """
