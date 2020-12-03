# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import tempfile
import typing as ty

import attr
from attr.validators import instance_of as is_a

import gffutils
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases import data

@attr.s()
class TranscriptInfo:
    so_rna_type = attr.ib(validator=is_a(str))
    regions: ty.List[data.SequenceRegion] = attr.ib(validator=is_a(list))


def load_coordinates(assembly_id: str, path: str) -> SqliteDict:
    mapping = SqliteDict()
    with tempfile.NamedTemporaryFile() as tmp:
        db = gffutils.create_db(path, tmp.name)
        for ncrna_gene in db.features_of_type('ncRNA_gene'):
            for transcript in db.children(ncrna_gene):
                exons: ty.List[data.Exon] = []
                gff_exons = db.children(transcript, feature_type='exon', order_by='start')
                for exon in exons:
                    exons.append(data.Exon(
                        start=exon.start,
                        stop=exon.stop
                    ))

                pid = transcript["ID"][0]
                if pid in mapping:
                    raise ValueError(f"Duplicate ensembl ID seen {pid}")
                mapping[pid] = TranscriptInfo(
                        so_rna_type=transcript.featuretype,
                        regions=[
                            data.SequenceRegion(
                                assembly_id=assembly_id,
                                chromosome=ncrna_gene.chrom,
                                strand=transcript.strand,
                                exons=exons,
                                coordinate_system=data.CoordinateSystem.one_based(),
                            )
                        ]
                    )
    return mapping
