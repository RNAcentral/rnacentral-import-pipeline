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
import nltk
import re


words = set(nltk.corpus.words.words())
ignore_ids = [
    "com", "mir-", "arg", "asn", "cys", "gln", "glu", "gly", "ile", "phe", "thr", "trp", "tyr", "val", "asx", "glx",
    "xaa"
]
words.update(ignore_ids)
special_char = re.compile('[@!#$%^&()<>?/\[\]\'}{~:]')
nts = re.compile('^[acgu]+$')


def check_id(item):
    if item.isnumeric() or item.lower() in words:
        result = None
    elif len(item) > 2 and not special_char.search(item) and not nts.search(item.lower()) and "\\" not in item:
        result = item
    else:
        result = None

    return result


@click.command()
@click.argument('database')
@click.argument('filename')
@click.argument('output')
def main(database, filename, output):
    """
    Check ids and create file that will be used by RNAcentral-references.
    """
    remove_dot = ["ensembl_gene", "ensembl_gencode_gene", "ensembl_metazoa_gene"]
    split_on_comma = ["flybase_gene_synonym", "pombase_gene_synonym", "refseq_gene_synonym"]

    with open(filename, 'r') as input_file:
        with open(output, 'w') as output_file:
            while line := input_file.readline():
                line = line.rstrip()
                line = line.split('|')

                if len(line) == 4:
                    get_gene = line[0]
                    get_primary_id = line[1]
                    urs = line[2]
                    taxid = line[3]

                    # remove "."
                    if database in remove_dot and "." in get_gene:
                        get_gene = get_gene.split('.')[0]

                    # split on ","
                    gene_results = []
                    if database in split_on_comma:
                        gene_list = get_gene.split(',')
                        for item in gene_list:
                            item = check_id(item)
                            if item:
                                gene_results.append(item)

                    if gene_results:
                        primary_id = check_id(get_primary_id)
                        for gene in gene_results:
                            if gene and primary_id and gene != primary_id:
                                output_file.write(gene + '|' + primary_id + '|' + urs + '_' + taxid + '\n')
                    else:
                        gene = check_id(get_gene)
                        primary_id = check_id(get_primary_id)
                        if gene and primary_id and gene != primary_id:
                            output_file.write(gene + '|' + primary_id + '|' + urs + '_' + taxid + '\n')

                else:
                    get_primary_id = line[0]
                    urs = line[1]
                    taxid = line[2]

                    # check if it is a valid id
                    primary_id = check_id(get_primary_id)

                    if primary_id:
                        output_file.write(primary_id + '|' + urs + '_' + taxid + '\n')


if __name__ == '__main__':
    main()
