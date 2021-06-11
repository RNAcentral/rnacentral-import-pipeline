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
import json
import itertools as it

import pymysql


DISALLOWED_DATABASE_TERMS = {
    "mirror",
    "chok1gs",  # skip golden hamster genome (same taxid as CriGri)
    "female",  # skip female naked mole rat genome
}


def database_key(name):
    parts = name.strip().split("_core_")
    numbers = parts[1].split("_")
    major = int(numbers[0])
    match = re.match(r"^(\d+)([^\d]?)$", numbers[1])
    if not match:
        raise ValueError("Cannot process: " + name)
    minor = int(match.group(1))
    suffix = ""
    if match.group(2):
        suffix = match.group(2)
    return (parts[0], major, minor, suffix)


def major(database):
    return database_key(database)[1]


def select_max(handle):
    databases = []
    for line in handle:
        if any(t in line for t in DISALLOWED_DATABASE_TERMS):
            continue
        if "mus_musculus" in line and line.count("_") != 4:
            continue  # skip mouse strains Mouse 129S1/SvImJ
        if "bacteria" in line:
            continue
        if "_core_" in line:
            databases.append(line.strip())

    max_major = max(database_key(d)[1] for d in databases)
    possible = filter(lambda d: major(d) == max_major, databases)
    grouped = it.groupby(possible, lambda d: database_key(d)[0])
    for _, databases in grouped:
        yield max(databases, key=database_key)


def databases(connection):
    cursor = connection.cursor()
    cursor.execute("show databases")
    db_names = [r[0] for r in cursor.fetchall()]
    cursor.close()
    databases = select_max(db_names)
    return databases


def run_queries(connection, query):
    for database in databases(connection):
        connection.select_db(database)
        cursor = connection.cursor(cursor=pymysql.cursors.DictCursor)
        cursor.execute(query)
        yield (database, cursor.fetchall())
        cursor.close()


def run_queries_across_databases(connection_handle, query_handle):
    specs = json.load(connection_handle)
    query = query_handle.read()
    for name, spec in specs.items():
        if not name.lower().startswith("ensembl"):
            continue
        if "command" in spec:
            del spec["command"]
        connection = pymysql.Connection(**spec)
        for results in run_queries(connection, query):
            yield results
