#!/usr/bin/env python
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
    This function creates a file to store URS and a file with job_ids|URS.
    These files can be used to create the metadata for the RNAcentral website.
    :param filename: file containing ids
    :param output: file to be created
    :return: None
    """
    type_id = output.split('_')[0]

    with open(filename, "r") as input_file:
        with open(output, 'w') as output_file:
            while line := input_file.readline():
                line = line.rstrip()
                line = line.split('|')

                if type_id == 'urs':
                    urs = line[-1]
                    output_file.write(urs + '\n')
                elif type_id == 'job' and len(line) == 2:
                    job = line[0]
                    urs = line[1]
                    if job and urs:
                        output_file.write(job + '|' + urs + '\n')
                elif type_id == 'job' and len(line) == 3:
                    job = line[0]
                    primary = line[1]
                    urs = line[2]
                    if job and urs:
                        output_file.write(job + '|' + urs + '\n')
                    if primary and urs:
                        output_file.write(primary + '|' + urs + '\n')


if __name__ == "__main__":
    main()
