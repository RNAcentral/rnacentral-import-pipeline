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

import typing as ty
from pathlib import Path

from rnacentral_pipeline.databases.data import TRnaScanResults


def parse(path: Path) -> ty.Iterator[TRnaScanResults]:
    """
    Parse a the results of tRNAScan-SE into an iterator of TRnaScanResults
    objects. This assumes the file is correclty formatted, ie has 3 header
    lines.
    """
    with path.open("r") as raw:
        _ = next(raw)
        _ = next(raw)
        sep = next(raw)
        assert sep.startswith("--")
        for line in raw:
            yield TRnaScanResults.from_line(line)
