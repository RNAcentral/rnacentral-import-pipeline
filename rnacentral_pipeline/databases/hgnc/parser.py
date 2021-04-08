# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import json
from pathlib import Path
import typing as ty

from rnacentral_pipeline.databases import data



def rnacentral_entry(entry: RawEntry, ensembl_mapping: EnsemblMapping) -> ty.Optional[data.Entry]
     """
     Map HGNC ncRNAs to RNAcentral using RefSeq, Vega, gtRNAdb accessions
     and sequence matches.
     """

     if entry.refseq_accession:
         rnacentral_id = self.get_rnacentral_id(entry['refseq_accession'][0])
         return rnacentral_id

     elif entry.gtrnadb_id:
         gtrnadb_id = raw_entry.gtrnadb_id
         if gtrnadb_id:
             rnacentral_id = gtrnadb_to_rnacentral(gtrnadb_id)

     elif entry.ensembl_gene_id:
         gene_id = entry.ensembl_gene_id
         fasta = ensembl_sequence(gene_id)
         if not fasta:
             return None

         md5 = get_md5(fasta)
         rnacentral_id = self.get_rnacentral_id_by_md5(md5)
         if rnacentral_id:
             return rnacentral_id
         else:
             return ensembl_mapping.get(gene_id, None)
    else:
        LOGGER.info("Cannot map %s", entry)
        return None


def parse(path: Path) -> ty.Iterable[data.Entry]:
    with path.open('r') as handle:
        raw_data = json.load(handle)

    for raw in raw_data['response']['docs']:
        raw_entry = RawEntry.from_raw(raw)
        mapped = rnacentral_entry(raw_entry)
        if mapped:
            yield mapped
