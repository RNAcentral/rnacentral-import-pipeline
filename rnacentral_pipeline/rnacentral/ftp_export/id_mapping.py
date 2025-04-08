# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import csv
import re
import typing as ty

from attr import frozen

from rnacentral_pipeline import psql


@frozen
class IdMapping:
    upi: str
    accession: str
    taxid: int
    external_id: str
    optional_id: str
    rna_type: str
    gene: str
    product: str
    database: str

    def entry(self) -> ty.List[str]:
        database = self.database
        accession = self.external_id
        gene = self.gene or ""
        gene = gene.replace("\t", " ")

        if self.database == "PDBE":
            database = "PDB"
            accession = "%s_%s" % (self.external_id, self.optional_id)
        elif self.database == "HGNC":
            accession = self.accession
        elif self.database == "ENA":
            if self.rna_type == "piRNA":
                gene = self.product
            accession = self.accession
        elif self.database == "MIRBASE":
            gene = self.optional_id
        elif self.database == "ENSEMBL":
            gene = self.optional_id
        elif self.database == "RFAM":
            if self.rna_type == "pre_miRNA":
                acc_range, endpoints, _ = self.accession.split(":", 2)
                acc = re.sub(r"\.\d\d+$", "", acc_range)
                start, stop = endpoints.split("..", 1)
                gene = f"{acc}/{start}-{stop}"

        return [
            self.upi,
            database,
            accession,
            str(self.taxid),
            self.rna_type,
            gene,
        ]


def generate_file(json_file, output):
    """
    This will generate a TSV mapping file given the input TSV.
    """

    writer = csv.writer(output, delimiter="\t")
    for entry in psql.json_handler(json_file):
        data = IdMapping(**entry)
        writer.writerow(data.entry())
