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

import click
import psycopg2
import psycopg2.extras
import pandas as pd


@click.group("scan-imports")
def cli():
    """
    A group of commands to scan imports and decide what to run
    """
    pass

"""
This is the process I think

manual selection -> csv file
csv file: db_name, remote -> nf runs process(check_db_md5) -> list of db_name: md5
list of db_name: md5 -> nf runs process(select_for_import) -> db_selection.config

main pipeline includes selection.config to switch on/off the right dbs

md5 creation can be done in shell for now, could do something with stripping the
metadata of json files to compare only the actual data md5 (date would change the overall sum)

"""




@cli.command("select-for-import")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("db_md5_map")
@click.argument("output", default="db_selection.config")
def select_db_to_import(db_md5_map, output, db_url=None, type=click.Path(writable=True,dir_okay=False,file_okay=True)):
    """
    Takes the map of db name to md5 sum and queries our DB to select those DBs that can usefully be imported

    Outputs a config file that the weekly import includes to switch on/off the relevant DBs
    """

    selection_template = """params {{
        databases {{
            {0}
            {1}
        }}
    }}"""

    latest_checksums = pd.read_csv(db_md5_map, names=["db_name", "checksum"])

    conn = psycopg2.connect(db_url)
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    cur.execute("SELECT * FROM rnc_import_tracker;")

    prev_checksums = pd.DataFrame(cur.fetchall())

    selection = latest_checksums.join(prev_checksums.set_index("db_name"), on="db_name").query("checksum != file_md5")['db_name'].values
    deselection = latest_checksums.join(prev_checksums.set_index("db_name"), on="db_name").query("checksum == file_md5")['db_name'].values

    selection = [f"{s}.run = true" for s in selection]
    deselection = [f"{s}.run = false" for s in deselection]

    activation_string = "\n\t\t".join(selection)
    deactivation_string = "\n\t\t".join(deselection)

    with open(output, 'w') as selection_config:
        selection_config.write(selection_template.format(activation_string, deactivation_string))
