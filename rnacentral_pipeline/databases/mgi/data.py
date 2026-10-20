# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from __future__ import annotations

import logging
import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline import db

LOGGER = logging.getLogger(__name__)


def split_ids(row, key: str) -> ty.List[str]:
    """
    MGI packs repeated ids into one '|' separated column.
    """
    value = row.get(key) or ""
    return [part for part in value.split("|") if part]


@attr.s(frozen=True)
class MgiEntry:
    """
    One marker from MRK_Sequence.rpt, reduced to what we can import.
    """

    mgi_id: str = attr.ib(validator=is_a(str))
    symbol: str = attr.ib(validator=is_a(str))
    name: str = attr.ib(validator=is_a(str))
    feature_type: str = attr.ib(validator=is_a(str))
    chromosome: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))
    refseq_ids: ty.Tuple[str, ...] = attr.ib(validator=is_a(tuple))
    ensembl_ids: ty.Tuple[str, ...] = attr.ib(validator=is_a(tuple))

    @classmethod
    def from_row(cls, row: ty.Dict[str, str]) -> MgiEntry:
        chromosome = row["chromosome"] or None
        if chromosome == "UN":
            chromosome = None

        return cls(
            mgi_id=row["mgi_marker_accession_id"],
            symbol=row["marker_symbol"],
            name=row["marker_name"],
            feature_type=row["feature_type"],
            chromosome=chromosome,
            refseq_ids=tuple(split_ids(row, "refseq_transcript_ids")),
            ensembl_ids=tuple(split_ids(row, "ensembl_transcript_ids")),
        )


@attr.s()
class Context:
    """
    Everything needed to map MGI markers to RNAcentral, resolved up front.

    MGI publishes gene models, not sequences, so an entry only exists for us if
    one of its transcripts is already in RNAcentral. All the lookups happen in
    bulk when this is built, so mapping an individual marker is a dict lookup.
    """

    refseq: ty.Dict[str, str] = attr.ib(factory=dict)
    ensembl: ty.Dict[str, str] = attr.ib(factory=dict)
    sequences: ty.Dict[str, str] = attr.ib(factory=dict)

    @classmethod
    def build(cls, db_url: str, entries: ty.List[MgiEntry]) -> Context:
        # Imported here as helpers needs MgiEntry from this module.
        from rnacentral_pipeline.databases.mgi import helpers

        conn = db.connect(db_url)
        try:
            context = cls(
                refseq=helpers.refseq_mapping(
                    conn, [i for e in entries for i in e.refseq_ids]
                ),
                ensembl=helpers.ensembl_mapping(
                    conn, [i for e in entries for i in e.ensembl_ids]
                ),
            )
            LOGGER.info(
                "Resolved %i refseq and %i ensembl transcripts",
                len(context.refseq),
                len(context.ensembl),
            )

            resolved = [context.urs_for(e) for e in entries]
            context.sequences = helpers.sequence_mapping(
                conn, [urs for urs in resolved if urs]
            )
        finally:
            conn.close()
        return context

    def urs_for(self, entry: MgiEntry) -> ty.Optional[str]:
        """
        Map a single marker to a URS, preferring RefSeq over Ensembl.
        """
        for refseq_id in entry.refseq_ids:
            urs = self.refseq.get(refseq_id)
            if urs:
                return urs

        for ensembl_id in entry.ensembl_ids:
            urs = self.ensembl.get(ensembl_id)
            if urs:
                return urs

        return None
