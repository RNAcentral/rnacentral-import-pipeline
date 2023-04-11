#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright [2009-present] EMBL-European Bioinformatics Institute
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


@click.command()
@click.argument('filename')
@click.argument('output')
def main(filename, output):
    """
    This function creates a file with the metadata of a given database.
    :param filename: file containing ids
    :param output: file to be created
    :return: None
    """
    with open(filename, "r") as input_file:
        database = filename.split(".")[0]

        with open(output, "w") as output_file:

            while line := input_file.readline():
                line = line.rstrip()
                line = line.split('|')

                if len(line) < 3:
                    job_id = line[0].lower()
                    urs = line[1].lower()

                    if job_id:
                        output_file.write(job_id + "|" + database + "\n")
                    if urs:
                        output_file.write(urs + "|" + "rnacentral" + "\n")
                    if job_id and urs:
                        output_file.write(job_id + "|" + "rnacentral" + "|" + urs + "\n")

                else:
                    job_id = line[0].lower()
                    primary_id = line[1].lower()
                    urs = line[2].lower()

                    if urs:
                        output_file.write(urs + "|" + "rnacentral" + "\n")
                    if primary_id:
                        output_file.write(primary_id + "|" + database + "\n")
                    if primary_id and urs:
                        output_file.write(primary_id + "|" + "rnacentral" + "|" + urs + "\n")
                    if job_id and urs:
                        output_file.write(job_id + "|" + "rnacentral" + "|" + urs + "\n")
                    if job_id and primary_id:
                        output_file.write(job_id + "|" + database + "|" + primary_id + "\n")


if __name__ == "__main__":
    main()
