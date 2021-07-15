# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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
from pathlib import Path
import typing as ty

from rnacentral_pipeline.cpat.data import CpatResult
from rnacentral_pipeline.cpat.data import CpatCutoffs


def cutoffs(directory: Path) -> CpatCutoffs:
    cutoffs = CpatCutoffs()
    for filename in directory.glob("dat/*_cutoff.txt"):
        source = filename.name.split("_", 1)[0].lower()
        with open(filename, 'r') as raw:
            lines = raw.readlines()
            assert len(lines) == 2, f"Unexpected cutoff format {filename}"
            if not lines[0].startswith("Coding Probability Cutoff: "):
                raise ValueError(f"Bad config format in {filename}")
            cutoff = float(lines[0].split(": ")[1])
            cutoffs.add_cutoff(source, cutoff)
    if not cutoffs.cutoffs:
        raise ValueError(f"Loaded no data from {directory}")
    return cutoffs


def parse(cutoffs: CpatCutoffs, model_name: str, results: Path) -> ty.Iterable[CpatResult]:
    with results.open('r') as raw:
        reader = csv.DictReader(raw, delimiter='\t')
        for data in reader:
            yield CpatResult.build(data, model_name, cutoffs)
