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

import re
import os
import csv
import logging
from glob import glob

import typing as ty
from pathlib import Path

from Bio import SeqIO

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from . import ribotyper
from .data import TravelerResult

LOGGER = logging.getLogger(__name__)


def ribotyper_models(directory, colored=True, allow_missing=False):
    """
    Look at the files in the given directory to find all the URS-model pairs
    that have been computed.
    """

    dir_path = Path(directory)
    ribo_results = ribotyper.as_dict(dir_path)
    seen = False
    for filename in dir_path.glob('*.fasta'):
        basename = os.path.basename(filename)
        pair, _ = os.path.splitext(basename)
        urs, model = pair.split('-', 1)
        ribo_result = ribo_results[urs]
        result = TravelerResult.build(
            urs,
            model,
            directory,
            ribo_result,
            colored=colored,
        )
        if result.is_valid():
            seen = True
            yield result

    if not seen:
        msg = "Found no possible models in: %s" % directory
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.error(msg)


def rfam_models(directory, colored=True, allow_missing=False):
    seen = False
    pattern = '*.svg'
    if colored:
        pattern = '*.colored.svg'
    for model_directory in glob(os.path.join(directory, '*')):
        model = os.path.basename(model_directory)
        for filename in glob(os.path.join(model_directory, pattern)):
            if not colored and '.colored' in filename:
                continue
            basename = os.path.basename(filename)
            urs, _ = os.path.splitext(basename)
            urs = urs.replace('.colored', '')
            result = TravelerResult.build(
                urs,
                model,
                directory,
                None,
                colored=colored,
                is_rfam=True,
            )
            if result.is_valid():
                seen = True
                yield result

    if not seen:
        msg = "Found no possible models in: %s" % directory
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.error(msg)


def models(directory, colored=True, allow_missing=False):
    has_fasta = bool(glob(os.path.join(directory, '*.fasta')))
    if has_fasta:
        return ribotyper_models(
            directory, 
            colored=colored,
            allow_missing=allow_missing,
        )

    has_rfam = bool(glob(os.path.join(directory, 'RF*')))
    if has_rfam:
        return rfam_models(
            directory, 
            colored=colored,
            allow_missing=allow_missing,
        )

    raise ValueError("Do not know how to parse contents of: %s" % directory)
