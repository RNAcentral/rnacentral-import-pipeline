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
from datetime import date

from lxml import etree
from lxml.builder import E

from .data import builder as raw_builder


def write_entries(handle, results):
    """
    This will write all entries in the given results iterable to the given
    handle. It then returns the number of entries written.
    """

    count = 0
    for result in results:
        count += 1
        handle.write(etree.tostring(result).decode())
        handle.write("\n")
    return count


def write(results, handle, count_handle):
    """
    This will create the required root XML element and place all the given
    XmlEntry objects as ElementTree.Element's in it. This then produces the
    string representation of that document which can be saved.
    """

    # pylint: disable=no-member
    handle.write("<database>")
    handle.write(etree.tostring(E.name("RNAcentral")).decode())
    handle.write(
        etree.tostring(
            E.description("a database for non-protein coding RNA sequences")
        ).decode()
    )
    handle.write(etree.tostring(E.release("1.0")).decode())
    handle.write(
        etree.tostring(E.release_date(date.today().strftime("%d/%m/%Y"))).decode()
    )

    handle.write("<entries>")
    count = write_entries(handle, results)
    handle.write("</entries>")

    if not count:
        raise ValueError("No entries found")

    count_handle.write(str(count))
    handle.write(etree.tostring(E.entry_count(str(count))).decode())
    handle.write("</database>")


def parse(input_handle):
    for line in input_handle:
        if not line:
            continue
        entry = json.loads(line)
        yield raw_builder(entry)


def as_xml(input_handle, output_handle, count_handle):
    data = parse(input_handle)
    write(data, output_handle, count_handle)


def release_note(output_handle, release, count_handles):
    total = 0
    for handle in count_handles:
        total += sum(int(line) for line in handle)
    now = date.today().strftime("%d-%B-%Y")
    output_handle.write("release=%s\n" % release)
    output_handle.write("release_date=%s\n" % now)
    output_handle.write("entries=%s\n" % total)
