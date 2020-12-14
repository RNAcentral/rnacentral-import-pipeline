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

import typing as ty

from rnacentral_pipeline.databases import data


def as_exon(location) -> data.Exon:
    """
    Build an Exon from a biopython location object.
    """

    return data.Exon(start=location.start + 1, stop=int(location.end))


def as_exons(feature) -> ty.Tuple[str, ty.List[data.Exon]]:
    """
    Determine all Exons in this feature.
    """

    parts = [feature.location]
    if hasattr(feature.location, 'parts'):
        parts = feature.location.parts

    strand = parts[0].strand
    return strand, [as_exon(l) for l in parts]


def regions(record, feature) -> ty.List[data.SequenceRegion]:
    accessions = record.annotations['accessions']
    if len(accessions) > 1 and accessions[0][-1] == '-':
        accessions = [''.join(accessions)]
    acc = accessions[0]
    assert ':' in acc, "Invalid accession (%s) for %s" % (acc, record)
    parts = acc.split(':')
    assembly_id = parts[1]
    chromosome_name = parts[2]

    strand, exons = as_exons(feature)
    return [data.SequenceRegion(
        chromosome=chromosome_name,
        strand=strand,
        exons=exons,
        assembly_id=assembly_id,
        coordinate_system=data.CoordinateSystem.one_based(),
    )]


