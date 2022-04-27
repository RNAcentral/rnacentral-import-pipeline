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


import os
import tempfile
import subprocess

from io import StringIO

from rnacentral_pipeline import psql


def run_with_buffer(path, *replacements):
    """
    Run a query from the given file with the given replacements and return a
    io.StringIO object of the output.
    """

    with tempfile.NamedTemporaryFile("w") as tmp:
        with open(path, "r") as raw:
            query = raw.read()
            for (initial, replacement) in replacements:
                if initial.startswith(":'"):
                    replacement = "'" + replacement + "'"
                query = query.replace(initial, replacement)
            tmp.write(query)
            tmp.flush()

        cmd = subprocess.run(
            ["psql", "-f", tmp.name, os.environ["PGDATABASE"]],
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        cmd.check_returncode()
        return StringIO(cmd.stdout)


def run_with_replacements(path, *replacements, take_all=False, **kwargs):
    """
    Run a query from the given file with the given replacements.
    """

    buf = run_with_buffer(path, *replacements)
    results = psql.json_handler(buf)

    try:
        if take_all:
            return list(results)
        return next(results)
    except StopIteration:
        raise ValueError("Found nothing for %s" % str(replacements))


def run_range_as_single(upi, path):
    """
    Run a query through psql and alter a rna.id BETWEEN statement with a upi
    and taxid paramenter.
    """

    parts = upi.split("_")
    upi, taxid = parts[0], int(parts[1])
    return run_with_replacements(
        path,
        (
            "rna.id BETWEEN :min AND :max",
            "xref.upi ='%s' AND xref.taxid = %i" % (upi, taxid),
        ),
    )


def run_with_upi_taxid_constraint(rna_id, path, **kwargs):
    upi, taxid = rna_id.split("_")
    taxid = int(taxid)
    return run_with_replacements(
        path,
        (
            "where\n",
            """where
            xref.upi = '%s' and xref.taxid = %i
            and"""
            % (upi, taxid),
        ),
        **kwargs
    )


def compare_dicts(dict1, dict2, ignore_keys=[]):
    """
    Compare two dictionaries with option to ignore some ignore_keys

    Use this when some things might change, e.g. when getting data from a URL

    Taken from
    https://stackoverflow.com/questions/10480806/compare-dictionaries-ignoring-specific-keys
    """
    print(dict1)
    dict1_filt = {k:v for k,v in dict1.items() if k not in ignore_keys}
    dict2_filt = {k:v for k,v in dict2.items() if k not in ignore_keys}
    return dict1_filt == dict2_filt
