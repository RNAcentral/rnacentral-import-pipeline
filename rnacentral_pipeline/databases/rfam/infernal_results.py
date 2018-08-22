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

import re
import csv

import attr

CSV_COLUMNS = [
    'seq_name',
    'seq_from',
    'seq_to',
    'strand',
    'rfam_acc',
    'mdl_from',
    'mdl_to',
    'overlap',
    'e_value',
    'score',
]


def as_0_based(raw):
    """
    Turn one based position into into 0 based coordinates.
    """
    return int(raw) - 1


def convert_strand(raw):
    """
    Turn a '+'/'-' strand to the -1/1 that we use internally.
    """

    if raw == '+' or raw == 1:
        return 1
    if raw == '-' or raw == -1:
        return -1
    raise Exception("Unknown strand %s" % raw)


def dash_as_none(handle=lambda r: r):
    """
    Turn a dash into None and if there is a value use handle to convert it's
    value.
    """

    def func(raw):
        """
        The actual converter function.
        """

        if raw is None or raw == '-':
            return None
        return handle(raw)
    return func


def convert_trunc(raw):
    """
    Convert the turnc value to a tuple of (True, False) to indicate if it is 5'
    or 3' truncated.
    """

    if isinstance(raw, tuple):
        return raw
    if raw == 'no':
        return (False, False)
    if raw == "5'":
        return (True, False)
    if raw == "3'":
        return (False, True)
    if raw == "5'&3'":
        return (True, True)
    raise Exception("Unknown yes/no %s" % raw)


def convert_overlap(raw):
    """
    Convert the short overlap strings from infernal into more readable names.
    """

    if raw == '*' or raw == 'unique':
        return u'unique'
    if raw == '^' or raw == 'best':
        return u'best'
    if raw == '=' or raw == 'secondary':
        return u'secondary'
    raise Exception("Unknown overlap symbol %s" % raw)


@attr.s(frozen=True, slots=True)
class RfamHit(object):
    """
    This represents the result of an Rfam scan using infernal without any clan
    competition.
    """

    target_name = attr.ib()
    seq_acc = attr.ib(convert=dash_as_none())
    rfam_name = attr.ib()
    rfam_acc = attr.ib()
    mdl = attr.ib()
    mdl_from = attr.ib(convert=as_0_based)
    mdl_to = attr.ib(convert=int)
    seq_from = attr.ib(convert=as_0_based)
    seq_to = attr.ib(convert=int)
    strand = attr.ib(convert=convert_strand)
    trunc = attr.ib(convert=convert_trunc)
    infernal_pass = attr.ib(convert=int)
    infernal_gc = attr.ib(convert=float)
    bias = attr.ib(convert=float)
    score = attr.ib(convert=float)
    e_value = attr.ib(convert=float)
    inc = attr.ib()
    description = attr.ib(convert=dash_as_none())


@attr.s(frozen=True, slots=True)  # pylint: disable=R0903
class RfamClanHit(object):
    """
    This represents the result of an rfam scan using infernal after clan
    competition.
    """

    idx = attr.ib(convert=as_0_based)
    rfam_name = attr.ib()
    rfam_acc = attr.ib()
    seq_name = attr.ib()
    seq_acc = attr.ib(convert=dash_as_none())
    clan_name = attr.ib(convert=dash_as_none())
    mdl = attr.ib()
    mdl_from = attr.ib(convert=as_0_based)
    mdl_to = attr.ib(convert=int)
    seq_from = attr.ib(convert=as_0_based)
    seq_to = attr.ib(convert=int)
    strand = attr.ib(convert=convert_strand)
    trunc = attr.ib(convert=convert_trunc)
    infernal_pass = attr.ib(convert=int)
    infernal_gc = attr.ib(convert=float)
    bias = attr.ib(convert=float)
    score = attr.ib(convert=float)
    e_value = attr.ib(convert=float)
    inc = attr.ib()
    overlap = attr.ib(convert=convert_overlap)
    anyidx = attr.ib(convert=dash_as_none(as_0_based))
    afrct1 = attr.ib(convert=dash_as_none(float))
    afrct2 = attr.ib(convert=dash_as_none(float))
    winidx = attr.ib(convert=dash_as_none())
    wfrct1 = attr.ib(convert=dash_as_none())
    wfrct2 = attr.ib(convert=dash_as_none())
    description = attr.ib(convert=dash_as_none())


def parse(raw, clan_competition=False):
    """
    This can parse the results of an infernal scan. It can handle both the
    standard tbl format as well as the tbl produced after clan competition. To
    handle clan copetition add the clan_competition=True argument. In this case
    the function will return an instance of the RfamClanHit class instead of
    the RfamHit class.
    """

    model = RfamHit
    if clan_competition:
        model = RfamClanHit
    fields = attr.fields(model)
    count = len(fields) - 1

    for line in raw:
        if line.startswith('#'):
            continue
        parts = re.split(r'\s+', line.strip(), count)
        yield model(*parts)  # pylint: disable=star-args


def as_csv(tblout, output):
    """
    This will parse the givn tblout filehandle and turn the data into a CSV
    file that can be loaded into the database.
    """

    writer = csv.DictWriter(
        output,
        CSV_COLUMNS,
        extrasaction='ignore',
        delimiter=',',
        quotechar='"',
        quoting=csv.QUOTE_ALL,
        lineterminator='\n',
    )

    for hit in parse(tblout, clan_competition=True):
        if hit.overlap == 'unique' or hit.overlap == 'best':
            writer.writerow(attr.asdict(hit))
