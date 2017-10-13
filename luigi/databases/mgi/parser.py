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

This contains the logic for parsing MGI data files and producing Entry objects
for export to usable flat files.
"""

import csv

from database.data import Entry

from .helpers import primary_id
from .helpers import accession
from .helpers import taxon_id
from .helpers import exon
from .helpers import infer_rna_type
from .helpers import xref_data
from .helpers import chromosome
from .helpers import lineage
from .helpers import gene
from .helpers import name
from .helpers import symbol
from .helpers import start
from .helpers import stop
from .helpers import references
from .helpers import species
from .helpers import common_name


def lines(raw):
    """
    Produces an iterable of all ines in the file. This will correct the issues
    with header being over 2 lines so a normal CSV parser can parse the file.
    """

    header = '\t'.join([next(raw).strip(), next(raw).strip()])
    header = header.lower()
    yield header.replace(' ', '_')
    for line in raw:
        yield line


def parser(filename):
    """
    Parses the file and produces an iterable of all entries as MGI objects.
    """

    with open(filename, 'rb') as raw:
        for data in csv.DictReader(lines(raw), delimiter='\t'):
            yield Entry(
                primary_id=primary_id(data),
                accession=accession(data),
                ncbi_tax_id=taxon_id(data),
                database='MGI',
                sequence='',
                exons=exon(data),
                rna_type=infer_rna_type(data) or '',
                url='',
                division='MUS',
                is_composite='N',
                xref_data=xref_data(data),
                chromosome=chromosome(data),
                species=species(data),
                common_name=common_name(data),
                lineage=lineage(data),
                gene=gene(data),
                optional_id=symbol(data),
                description=name(data),
                seq_version='1',
                feature_location_start=start(data),
                feature_location_end=stop(data),
                references=references(data),
            )


def rna_entries(filename):
    """
    Parses the file and produces an iterable of only the RNA entries.
    """

    for entry in parser(filename):
        if entry.rna_type:
            yield entry
