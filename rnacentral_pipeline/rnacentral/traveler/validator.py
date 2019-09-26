# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import six

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

import numpy as np
from scipy import stats

from rnacentral_pipeline import psql


REJECTED_RRNAS = {
    'group_I_intron',
    'group_II_intron',
    'large_subunit_rRNA',
}

ALLOWED_RFAM = {
    'RNase_MRP_RNA',
    'RNase_P_RNA',
    'SRP_RNA',
    'Y_RNA',
    'antisense_RNA',
    'autocatalytically_spliced_intron',
    'hammerhead_ribozyme',
    'ncRNA',
    'pre_miRNA',
    'rRNA',
    'ribozyme',
    'snRNA',
    'snoRNA',
    'tRNA',
    'telomerase_RNA',
    'tmRNA',
    'vault_RNA',
}


@attr.s()
class ShouldShow(object):
    urs = attr.ib(validator=is_a(six.text_type))
    model_type = attr.ib(validator=is_a(six.text_type))
    rna_type = attr.ib(validator=is_a(six.text_type))
    overlap_count = attr.ib(validator=is_a(int))
    zscore = attr.ib(validator=optional(is_a(float)))

    @classmethod
    def build(cls, data, zscore):
        return cls(
            urs=data['urs'],
            model_type=data['model_type'],
            rna_type=data['rna_type'],
            overlap_count=int(data['overlap_count']),
            zscore=zscore,
        )

    @property
    def should_show(self):
        if self.model_type == 'rfam':
            return self.__rfam_should_show__()
        if self.model_type == 'rrna':
            return self.__rrna_should_show__()
        raise ValueError("Unknown model type %s" % self.model_type)

    def __rfam_should_show__(self):
        if self.rna_type == 'lnc_RNA':
            return False
        if self.rna_type in ALLOWED_RFAM:
            return self.zscore < 10
        raise ValueError("Unknown Rfam model %s" % self)

    def __rrna_should_show__(self):
        if self.rna_type == 'small_subunit_rRNA':
            return self.overlap_count < 14
        if self.rna_type == 'rRNA_5S':
            return self.overlap_count < 3
        if self.rna_type in REJECTED_RRNAS:
            return False
        raise ValueError("Unknown type of rRNA %s" % self)

    def writeable(self):
        return [
            self.urs,
            self.zscore,
            self.should_show,
        ]


def parse_results(raw):
    return list(psql.json_handler(raw))


def should_show(data):
    """
    Compute the should_show values for all entries int he given file assuming
    the given RNA type.
    """

    sequence_lengths = np.array([d['sequence_length'] for d in data], dtype=np.float)
    model_lengths = np.array([d['model_length'] for d in data], dtype=np.float)

    zscores = np.array([np.nan for i in range(len(sequence_lengths))])
    if len(model_lengths) == len(sequence_lengths):
        zscores = stats.zscore(sequence_lengths / model_lengths)

    for index, entry in enumerate(data):
        score = zscores[index]
        if not np.isnan(score):
            score = None
        yield ShouldShow.build(entry, score)


def from_file(raw):
    return should_show(parse_results(raw))
