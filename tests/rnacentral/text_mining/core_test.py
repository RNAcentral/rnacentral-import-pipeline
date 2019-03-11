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

import six
import attr
import pytest
import textblob as tb

from rnacentral_pipeline.rnacentral.text_mining import core


def test_can_find_extact_names_in_text():
    matcher = core.NameMatcher.build(['HOTAIR', 'XIST'])
    selector = core.SentenceSelector.build()
    with open('data/text-mining/PMID:29141248.txt', 'r') as raw:
        blob = tb.TextBlob(raw.read())
        assert [attr.asdict(m) for m in core.matches(blob, selector, matcher)] == [
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'BACKGROUND/AIMS:The HOX transcript antisense intergenic RNA (HOTAIR), a long '
                    'non-coding RNA (lncRNA), plays an important role in the '
                    'pathogenesis and progression of multiple tumors.'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'The aim of the present study was to evaluate whether common '
                    'single nucleotide polymorphisms (SNPs) in HOTAIR are related '
                    'to hepatocellular carcinoma (HCC) susceptibility in a '
                    'Chinese population.'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'METHODS:We genotyped three SNPs of HOTAIR in a hepatocellular '
                    'carcinoma (HCC) case-control study, including 482 cases and '
                    '520 control subjects.'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'The allele-specific effects on HOTAIR expression in HCC were '
                    'confirmed by real time quantitative PCR and luciferase '
                    'activity assays.'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'The influence of HOTAIR SNPs on the proliferation of HCC '
                    'cells was evaluated using a CCK-8 assay.'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'RESULTS:Significant associations were observed between the HOTAIR '
                    'rs920778 C>T polymorphism and HCC risk (TT versus CC: OR = '
                    '1.634, 95% CI =1.028-2.598, P = 0.046) and the '
                    'allelic model (allele T versus allele C: OR =1.293, 95% CI = '
                    '1.060-1.577, P = 0.011).'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'RT-PCR and luciferase activity assay confirmed that '
                    'the rs920778 TT genotype induced significantly higher HOTAIR levels than did the '
                    'CC genotype (P < 0.05).'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'CONCLUSION:These results suggest that SNP rs920778 of HOTAIR acts as a '
                    'potential biomarker for predicting hepatocellular carcinoma, '
                    'and further studies are warranted to confirm these findings.'
                ),
                matches=[core.WordMatch('HOTAIR', 'HOTAIR')]
            )),
        ]


def test_can_find_expected_patterns_in_text():
    matcher = core.PatternMatcher.build([r'(?P<gene>HO[TX])', r'(?P<allele>rs[41]\d+)'])
    selector = core.SentenceSelector.build()
    with open('data/text-mining/PMID:29141248.txt', 'r') as raw:
        blob = tb.TextBlob(raw.read())
        assert [attr.asdict(m) for m in core.matches(blob, selector, matcher)] == [
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'BACKGROUND/AIMS:The HOX transcript antisense intergenic RNA (HOTAIR), a long '
                    'non-coding RNA (lncRNA), plays an important role in the '
                    'pathogenesis and progression of multiple tumors.'
                ),
                matches=[core.WordMatch('HOX', 'gene')]
            )),
            attr.asdict(core.MatchingSentence(
                sentence=tb.Sentence(
                    'However, no statistically significant differences of '
                    'rs4759314 and rs1899663 genotypes were observed between '
                    'patients and controls (both P > 0.05).'
                ),
                matches=[
                    core.WordMatch('rs4759314', 'allele'),
                    core.WordMatch('rs1899663', 'allele'),
                ]
            )),
        ]
