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
    This function creates a file for each database containing the manually annotated references
    :param filename: file containing ids
    :param output: file to be created
    :return: None
    """
    name = output.split("*")[0]

    with open(filename, "r") as input_file:
        with open(name + "hgnc", 'w') as hgnc, open(name + "pombase", 'w') as pombase, open(name + "sgd", 'w') as sgd, \
                open(name + "tair", 'w') as tair, open(name + "zfin", 'w') as zfin:
            while line := input_file.readline():
                line = line.rstrip()
                line = line.split('|')
                urs = line[0]
                database = line[1]
                pmid = line[2]
                doi = line[3]
                pmcid = line[4]

                if database.lower() == "hgnc":
                    hgnc.write(urs + '|' + pmid + '|' + doi + '|' + pmcid + '\n')
                elif database.lower() == "pombase":
                    pombase.write(urs + '|' + pmid + '|' + doi + '|' + pmcid + '\n')
                elif database.lower() == "sgd":
                    sgd.write(urs + '|' + pmid + '|' + doi + '|' + pmcid + '\n')
                elif database.lower() == "tair":
                    tair.write(urs + '|' + pmid + '|' + doi + '|' + pmcid + '\n')
                elif database.lower() == "zfin":
                    zfin.write(urs + '|' + pmid + '|' + doi + '|' + pmcid + '\n')


if __name__ == "__main__":
    main()
