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

import json
import itertools as it
from datetime import date

from lxml import etree
from lxml.builder import E

from .data import builder


def write_entries(handle, results):
    count = 0
    for result in results:
        count += 1
        handle.write(etree.tostring(result))
    return count


def write(results, handle, count_handle):
    """
    This will create the required root XML element and place all the given
    XmlEntry objects as ElementTree.Element's in it. This then produces the
    string representation of that document which can be saved.
    """

    handle.write('<database>')
    handle.write(etree.tostring(E.name('RNAcentral')))
    handle.write(etree.tostring(E.description('a database for non-protein coding RNA sequences')))
    handle.write(etree.tostring(E.release('1.0')))
    handle.write(etree.tostring(E.release_date(date.today().strftime('%d/%m/%Y'))))

    handle.write('<entries>')
    count = write_entries(handle, results)
    handle.write('</entries>')

    if not count:
        raise ValueError("No entries found")
    count_handle.write(str(count))

    handle.write(etree.tostring(E.entry_count(str(count))))
    handle.write('</database>')


def as_xml(input_handle, output_handle, count_handle):
    def parser():
        for line in input_handle:
            line = line.replace('\\\\', '\\')
            yield json.loads(line)
    data = it.imap(builder, parser())
    write(data, output_handle, count_handle)


def release_note(output_handle, count_handles):
    total = 0
    for handle in count_handles:
        total += sum(int(line) for line in handle)
    now = date.today().strftime('%d-%B-%Y')
    output_handle.write("release=9.0\n")
    output_handle.write("release_date=%s\n" % now)
    output_handle.write("entries=%s\n" % total)
