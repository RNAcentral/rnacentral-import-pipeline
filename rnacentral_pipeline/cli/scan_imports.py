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
import pandas as pd
import psycopg2
import psycopg2.extras


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
def select_db_to_import(
    db_md5_map,
    output,
    db_url=None,
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
):
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

    cur.close()
    conn.close()

    selection = (
        latest_checksums.join(prev_checksums.set_index("db_name"), on="db_name")
        .query("checksum != file_md5")["db_name"]
        .values
    )
    deselection = (
        latest_checksums.join(prev_checksums.set_index("db_name"), on="db_name")
        .query("checksum == file_md5")["db_name"]
        .values
    )

    selection = [f"{s}.run = true" for s in selection]
    deselection = [f"{s}.run = false" for s in deselection]

    activation_string = "\n\t\t".join(selection)
    deactivation_string = "\n\t\t".join(deselection)

    with open(output, "w") as selection_config:
        selection_config.write(
            selection_template.format(activation_string, deactivation_string)
        )


@cli.command("update-tracker")
@click.argument("latest_md5s")
@click.option("--db-url", envvar="PGDATABASE")
def update_tracker(latest_md5s, db_url):
    """
    Record the checksum of each database's remote file in rnc_import_tracker.

    Every database in the file is written, not only the ones that changed: the tracker
    says what we last imported, so dropping the unchanged rows would make every
    database look new again on the next run.
    """
    latest = pd.read_csv(latest_md5s, names=["db_name", "checksum"])

    conn = psycopg2.connect(db_url)
    try:
        with conn.cursor() as cur:
            names = sorted({str(name).upper() for name in latest["db_name"]})
            cur.execute(
                "SELECT descr, id FROM rnc_database WHERE descr = ANY(%s)", (names,)
            )
            ids = dict(cur.fetchall())

            rows = [
                (str(name).lower(), ids[str(name).upper()], str(checksum))
                for name, checksum in zip(latest["db_name"], latest["checksum"])
                if str(name).upper() in ids
            ]
            if not rows:
                return

            # db_name carries no unique constraint, so replace this run's rows rather
            # than upserting -- and only this run's. Truncating the table would forget
            # every database the run did not cover.
            cur.execute(
                "DELETE FROM rnc_import_tracker WHERE db_name = ANY(%s)",
                ([row[0] for row in rows],),
            )
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO rnc_import_tracker
                    (db_name, db_id, last_import_date, file_md5)
                VALUES %s
                """,
                rows,
                template="(%s, %s, CURRENT_TIMESTAMP, %s)",
            )
        conn.commit()
    finally:
        conn.close()
