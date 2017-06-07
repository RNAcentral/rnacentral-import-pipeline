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

import logging
import operator as op

import attr
from Bio import SeqIO
import luigi

LOGGER = logging.getLogger(__name__)


def slicer(*keys):
    getter = op.itemgetter(*keys)
    def fn(data):
        return dict(zip(keys, getter(data)))
    return fn


class RfamHitsImporter(luigi.Task):
    def stored_id(self, model, data):
        query = model.objects.filter(**data)
        if not query.exists():
            model.create(**data)
        return query.get().id

    def parse_hits(self):
        return []

    def data(self):
        # Split file by Rfam model:
        # awk '{ print >> $1".txt" }' /hps/nobackup/production/xfam/rfam/rfam13_rnac/rnac_all.txt

        model_slice = slicer('rfam', 'mdl_from', 'mdl_to')
        seq_slice = slicer('urs', 'seq_from', 'seq_to')

        for hit in self.parse_hits():
            mdl_range = model_slice(hit)
            urs_range = seq_slice(hit)
            mdl_id = self.stored_id(RfamHitRegion, mdl_range)
            urs_id = self.stored_id(URS, seq_range)
            yield RfamModelHit(
                urs=hit['urs'],
                rfam_model_range_id=mdl_id,
                sequence_hit_range_id=urs_id,
                e_value=hit['e_value'],
                score=hit['score'],
            )

    def run(self):
        for entry in self.data():
            entry.save()
