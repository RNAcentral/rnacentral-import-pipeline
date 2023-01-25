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
import csv


@click.command()
@click.argument('filename')
@click.argument('output')
def main(filename, output):
    """
    The original CSV file has a column containing a list of PMIDs.
    This function will create a file where each line will have an
    organism and a pmid. This file will be used to load the data
    into the load_organism table.
    :param filename: file downloaded from https://organisms.jensenlab.org/
    :param output: file to be created
    :return: None
    """
    csv.field_size_limit(500000000)  # 500 megabytes

    with open(filename, "r") as input_file:
        with open(output, "w") as output_file:
            csv_reader = csv.reader(input_file, delimiter="\t")
            csv_writer = csv.writer(output_file)

            for item in csv_reader:
                organism = item[0]
                pmid_list = list(map(int, item[1].split(" ")))
                for pmid in pmid_list:
                    csv_writer.writerow([organism, pmid])


if __name__ == "__main__":
    main()
