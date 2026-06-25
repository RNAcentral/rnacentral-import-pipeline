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

import json
from pathlib import Path

from rnacentral_pipeline.repeatmasker import parser

# Two query sequences scanned; only the first has repeat hits. The second must
# still produce a presence row (has_repeats = False) for idempotency.
RM_OUT = """\
   SW   perc perc perc  query    position in query           matching  repeat        position in repeat
score   div. del. ins.  sequence  begin  end (left)   repeat        class/family   begin  end (left)  ID

  463   12.3  0.6  1.7  URS0000000001_9606    11   210 (90) +  AluY      SINE/Alu          1   200 (12)   1
  250    8.0  0.0  0.0  URS0000000001_9606   300   350 (0)  C  L1MA      LINE/L1         (5)  120    1     2
"""

FASTA = """\
>URS0000000001_9606
%s
>URS0000000002_9606
%s
""" % (
    "A" * 350,
    "C" * 100,
)


def _run(tmp_path: Path):
    out = tmp_path / "seq.fasta.out"
    out.write_text(RM_OUT)
    fasta = tmp_path / "seq.fasta"
    fasta.write_text(FASTA)
    return list(parser.parse(fasta, out))


def test_parses_hits_and_absence(tmp_path: Path):
    data = _run(tmp_path)
    by_urs = {result.urs_taxid: (result, regions) for result, regions in data}
    assert set(by_urs) == {"URS0000000001_9606", "URS0000000002_9606"}

    hit_result, hit_regions = by_urs["URS0000000001_9606"]
    assert hit_result.has_repeats is True
    assert hit_result.repeat_count == 2
    # covered = (210-10) + (350-299) = 200 + 51 = 251 of 350
    assert hit_result.repeat_coverage == round(251 / 350, 4)

    empty_result, empty_regions = by_urs["URS0000000002_9606"]
    assert empty_result.has_repeats is False
    assert empty_result.repeat_count == 0
    assert empty_regions == []


def test_feature_coordinates_and_strand(tmp_path: Path):
    data = _run(tmp_path)
    regions = next(
        regions for result, regions in data
        if result.urs_taxid == "URS0000000001_9606"
    )
    alu, line = regions
    assert (alu.start, alu.stop, alu.strand) == (10, 210, "+")
    assert (line.start, line.stop, line.strand) == (299, 350, "-")

    urs, taxid, start, stop, metadata = alu.writeable()
    assert (urs, taxid, start, stop) == ("URS0000000001", "9606", "10", "210")
    assert json.loads(metadata)["class"] == "SINE/Alu"
