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

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

SO_PATTERN = re.compile(r"^SO:\d+$")

INSDC_SO_MAPPING = {
    "RNase_MRP_RNA": "SO:0000385",
    "RNase_P_RNA": "SO:0000386",
    "SRP_RNA": "SO:0000590",
    "Y_RNA": "SO:0000405",
    "antisense_RNA": "SO:0000644",
    "autocatalytically_spliced_intron": "SO:0000588",
    "circRNA": "SO:0002291",
    "guide_RNA": "SO:0000602",
    "hammerhead_ribozyme": "SO:0000380",
    "lncRNA": "SO:0001877",
    "miRNA": "SO:0000276",
    "ncRNA": "SO:0000655",
    "misc_RNA": "SO:0000673",
    "other": "SO:0000655",
    "precursor_RNA": "SO:0000185",
    "piRNA": "SO:0001035",
    "rasiRNA": "SO:0000454",
    "ribozyme": "SO:0000374",
    "scRNA": "SO:0000013",
    "scaRNA": "SO:0002095",
    "sgRNA": "SO:0001998",
    "siRNA": "SO:0000646",
    "snRNA": "SO:0000274",
    "snoRNA": "SO:0000275",
    "telomerase_RNA": "SO:0000390",
    "tmRNA": "SO:0000584",
    "vault_RNA": "SO:0000404",
    "rRNA": "SO:0000252",
    "tRNA": "SO:0000253",
    "bidirectional_promoter_lncrna": "SO:0002185",
    "3prime_overlapping_ncrna": "SO:0002120",
    "pre_miRNA": "SO:0001244",
    "sRNA": "SO:0000655",
}

NORMALIZE_TO_INSDC = {
    "sRNA": "other",
    "bidirectional_promoter_lncrna": "lncRNA",
    "3prime_overlapping_ncrna": "other",
}

SO_INSDC_MAPPING = {v: k for k, v in INSDC_SO_MAPPING.items()}
SO_INSDC_MAPPING["SO:0000035"] = "ncRNA"
SO_INSDC_MAPPING["SO:0000077"] = "antisense_RNA"
SO_INSDC_MAPPING["SO:0000122"] = "other"
SO_INSDC_MAPPING["SO:0000204"] = "ncRNA"
SO_INSDC_MAPPING["SO:0000209"] = "precursor_RNA"
SO_INSDC_MAPPING["SO:0000233"] = "other"
SO_INSDC_MAPPING["SO:0000243"] = "other"
SO_INSDC_MAPPING["SO:0000370"] = "ncRNA"
SO_INSDC_MAPPING["SO:0000375"] = "rRNA"
SO_INSDC_MAPPING["SO:0000376"] = "other"
SO_INSDC_MAPPING["SO:0000377"] = "other"
SO_INSDC_MAPPING["SO:0000378"] = "other"
SO_INSDC_MAPPING["SO:0000379"] = "other"
SO_INSDC_MAPPING["SO:0000383"] = "other"
SO_INSDC_MAPPING["SO:0000384"] = "other"
SO_INSDC_MAPPING["SO:0000387"] = "other"
SO_INSDC_MAPPING["SO:0000389"] = "other"
SO_INSDC_MAPPING["SO:0000391"] = "snRNA"
SO_INSDC_MAPPING["SO:0000392"] = "snRNA"
SO_INSDC_MAPPING["SO:0000393"] = "snRNA"
SO_INSDC_MAPPING["SO:0000395"] = "snRNA"
SO_INSDC_MAPPING["SO:0000396"] = "snRNA"
SO_INSDC_MAPPING["SO:0000398"] = "snRNA"
SO_INSDC_MAPPING["SO:0000399"] = "snRNA"
SO_INSDC_MAPPING["SO:0000407"] = "rRNA"
SO_INSDC_MAPPING["SO:0000587"] = "ncRNA"
SO_INSDC_MAPPING["SO:0000593"] = "snoRNA"
SO_INSDC_MAPPING["SO:0000594"] = "snoRNA"
SO_INSDC_MAPPING["SO:0000603"] = "ncRNA"
SO_INSDC_MAPPING["SO:0000650"] = "rRNA"
SO_INSDC_MAPPING["SO:0000651"] = "rRNA"
SO_INSDC_MAPPING["SO:0000652"] = "rRNA"
SO_INSDC_MAPPING["SO:0000990"] = "other"
SO_INSDC_MAPPING["SO:0001000"] = "rRNA"
SO_INSDC_MAPPING["SO:0001001"] = "rRNA"
SO_INSDC_MAPPING["SO:0001032"] = "other"
SO_INSDC_MAPPING["SO:0001171"] = "rRNA"
SO_INSDC_MAPPING["SO:0001179"] = "snoRNA"
SO_INSDC_MAPPING["SO:0001263"] = "ncRNA"
SO_INSDC_MAPPING["SO:0001267"] = "snoRNA"
SO_INSDC_MAPPING["SO:0001268"] = "snRNA"
SO_INSDC_MAPPING["SO:0001272"] = "tRNA"
SO_INSDC_MAPPING["SO:0001459"] = "other"
SO_INSDC_MAPPING["SO:0001463"] = "lncRNA"
SO_INSDC_MAPPING["SO:0001637"] = "rRNA"
SO_INSDC_MAPPING["SO:0001643"] = "telomerase_RNA"
SO_INSDC_MAPPING["SO:0001904"] = "lncRNA"
SO_INSDC_MAPPING["SO:0002128"] = "rRNA"
SO_INSDC_MAPPING["SO:0002168"] = "other"
SO_INSDC_MAPPING["SO:0005836"] = "ncRNA"
SO_INSDC_MAPPING["SO:0005857"] = "tRNA"
SO_INSDC_MAPPING["SO:0000653"] = "rRNA"
SO_INSDC_MAPPING["SO:0000254"] = "tRNA"
SO_INSDC_MAPPING["SO:0001036"] = "tRNA"
SO_INSDC_MAPPING["SO:0000256"] = "tRNA"
SO_INSDC_MAPPING["SO:0000257"] = "tRNA"
SO_INSDC_MAPPING["SO:0000258"] = "tRNA"
SO_INSDC_MAPPING["SO:0000259"] = "tRNA"
SO_INSDC_MAPPING["SO:0000260"] = "tRNA"
SO_INSDC_MAPPING["SO:0000261"] = "tRNA"
SO_INSDC_MAPPING["SO:0000262"] = "tRNA"
SO_INSDC_MAPPING["SO:0000263"] = "tRNA"
SO_INSDC_MAPPING["SO:0000264"] = "tRNA"
SO_INSDC_MAPPING["SO:0000265"] = "tRNA"
SO_INSDC_MAPPING["SO:0000268"] = "tRNA"
SO_INSDC_MAPPING["SO:0000266"] = "tRNA"
SO_INSDC_MAPPING["SO:0000270"] = "tRNA"
SO_INSDC_MAPPING["SO:0000272"] = "tRNA"
SO_INSDC_MAPPING["SO:0000269"] = "tRNA"
SO_INSDC_MAPPING["SO:0000271"] = "tRNA"
SO_INSDC_MAPPING["SO:0000766"] = "tRNA"
SO_INSDC_MAPPING["SO:0001172"] = "tRNA"
SO_INSDC_MAPPING["SO:0002129"] = "tRNA"
SO_INSDC_MAPPING["SO:0000267"] = "tRNA"
SO_INSDC_MAPPING["SO:0000273"] = "tRNA"
SO_INSDC_MAPPING["SO:0002345"] = "rRNA"
SO_INSDC_MAPPING["SO:0002247"] = "other"


class UnxpectedRnaType(Exception):
    """
    Raised when the RNA type is not an SO term and cannot be converted to one.
    """

    pass


def optionally(instance_type, **kwargs):
    """
    Return an attribute that is either none or of the given type.
    """
    return attr.ib(validator=optional(is_a(instance_type)), default=None, **kwargs)


def possibly_empty(instance_type, **kwargs):
    """
    Return an attribute that defaults to being empty and must be of the given
    type.
    """

    factory = instance_type
    if hasattr(instance_type, "empty"):
        factory = instance_type.empty

    return attr.ib(
        validator=is_a(instance_type), default=attr.Factory(factory), **kwargs
    )


def matches_pattern(pattern):
    def fn(instance, attribute, value):
        if not re.match(pattern, value):
            raise TypeError(
                "Bad value (%s) for %s in %s" % (value, attribute, instance)
            )

    return fn


def as_so_term(rna_type: str) -> str:
    if re.match(SO_PATTERN, rna_type):
        return rna_type

    if rna_type not in INSDC_SO_MAPPING:
        raise UnxpectedRnaType(rna_type)
    return INSDC_SO_MAPPING[rna_type]
