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
import nltk
import re


nltk.download('words')
words = set(nltk.corpus.words.words())
ignore_ids = [
    "com", "mir-", "arg", "asn", "cys", "gln", "glu", "gly", "ile", "phe", "thr", "trp", "tyr", "val", "asx", "glx",
    "xaa"
]
words.update(ignore_ids)
special_char = re.compile('[@!#$%^&()<>?/\[\]\'}{~:]')
nts = re.compile('^[acgu]+$')
numbers_and_dash = re.compile('^\d+[\-]\d+$')  # do not use ids like 6-1, 260-1, etc


def check_id(item):
    if item.isnumeric() or item.lower() in words or numbers_and_dash.search(item):
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
    remove_dot = ["ensembl", "ensembl_gencode", "ensembl_metazoa"]
    split_on_comma = ["flybase", "hgnc", "pombase", "refseq"]
    rfam_ignore = [
        "30_255", "30_292", "5S_rRNA", "5_8S_rRNA", "6A", "6S", "7SK", "C4", "CRISPR-DR10", "CRISPR-DR11",
        "CRISPR-DR12", "CRISPR-DR13", "CRISPR-DR14", "CRISPR-DR15", "CRISPR-DR16", "CRISPR-DR17", "CRISPR-DR18",
        "CRISPR-DR19", "CRISPR-DR2", "CRISPR-DR20", "CRISPR-DR21", "CRISPR-DR22", "CRISPR-DR23", "CRISPR-DR24",
        "CRISPR-DR25", "CRISPR-DR26", "CRISPR-DR27", "CRISPR-DR28", "CRISPR-DR29", "CRISPR-DR3", "CRISPR-DR30",
        "CRISPR-DR31", "CRISPR-DR32", "CRISPR-DR33", "CRISPR-DR34", "CRISPR-DR35", "CRISPR-DR36", "CRISPR-DR37",
        "CRISPR-DR38", "CRISPR-DR39", "CRISPR-DR4", "CRISPR-DR40", "CRISPR-DR41", "CRISPR-DR42", "CRISPR-DR43",
        "CRISPR-DR44", "CRISPR-DR45", "CRISPR-DR46", "CRISPR-DR47", "CRISPR-DR48", "CRISPR-DR49", "CRISPR-DR5",
        "CRISPR-DR50", "CRISPR-DR51", "CRISPR-DR52", "CRISPR-DR53", "CRISPR-DR54", "CRISPR-DR55", "CRISPR-DR56",
        "CRISPR-DR57", "CRISPR-DR58", "CRISPR-DR6", "CRISPR-DR60", "CRISPR-DR61", "CRISPR-DR62", "CRISPR-DR63",
        "CRISPR-DR64", "CRISPR-DR65", "CRISPR-DR66", "CRISPR-DR7", "CRISPR-DR8", "CRISPR-DR9", "F6", "Hairpin",
        "Hairpin-meta1", "Hairpin-meta2", "Hatchet", "P1", "P10", "P11", "P13", "P14", "P15", "P17", "P18", "P2", "P24",
        "P26", "P27", "P31", "P33", "P34", "P35", "P36", "P37", "P4", "P5", "P6", "P8", "P9", "ROSE", "S35", "S414",
        "S774", "S808", "SAM", "SL1", "SL2", "U1", "U11", "U12", "U1_yeast", "U2", "U3", "U4", "U4atac", "U5", "U54",
        "U6", "U6atac", "U7", "U8", "VA", "csRNA", "drum", "g2", "pRNA", "sar", "sul1", "t44", "tRNA", "tRNA-Sec",
        "tmRNA", "tp2", "tracrRNA"
    ]

    with open(filename, 'r') as input_file:
        with open(output, 'w') as output_file:
            while line := input_file.readline():
                line = line.rstrip()
                line = line.split('|')
                urs = line[0]
                taxid = line[1]
                primary_id = check_id(line[2])
                if primary_id and database in remove_dot and "." in primary_id:
                    primary_id = primary_id.split('.')[0]

                if primary_id and line[3:]:
                    for item in line[3:]:
                        if item:
                            get_id = item
                        else:
                            continue

                        # ignore some optional_id from Rfam
                        if database == "rfam" and get_id in rfam_ignore:
                            output_file.write('|' + primary_id + '|' + urs + '_' + taxid + '\n')
                            continue

                        # remove "."
                        if database in remove_dot and "." in get_id:
                            get_id = get_id.split('.')[0]

                        # split on ","
                        results = []
                        if database in split_on_comma:
                            list_of_ids = get_id.split(',')
                            for elem in list_of_ids:
                                elem = check_id(elem)
                                if elem:
                                    results.append(elem)

                        if results:
                            for db_id in results:
                                if db_id != primary_id:
                                    output_file.write(db_id + '|' + primary_id + '|' + urs + '_' + taxid + '\n')
                        else:
                            db_id = check_id(get_id)
                            if db_id and db_id != primary_id:
                                output_file.write(db_id + '|' + primary_id + '|' + urs + '_' + taxid + '\n')
                elif primary_id:
                    output_file.write(primary_id + '|' + urs + '_' + taxid + '\n')


if __name__ == '__main__':
    main()
