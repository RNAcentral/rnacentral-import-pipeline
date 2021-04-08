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

import attr

# {
#   10   │         "date_approved_reserved": "2009-07-20",
#   11   │         "vega_id": "OTTHUMG00000183508",
#   12   │         "locus_group": "non-coding RNA",
#   13   │         "status": "Approved",
#   14   │         "alias_symbol": [
#   15   │           "FLJ23569"
#   16   │         ],
#   17   │         "_version_": 1696435635999473700,
#   18   │         "uuid": "ccf6f769-0fcc-467d-8c39-299191bdeacc",
#   19   │         "rna_central_id": [
#   20   │           "URS00007E4F6E"
#   21   │         ],
#   22   │         "prev_name": [
#   23   │           "non-protein coding RNA 181",
#   24   │           "A1BG antisense RNA (non-protein coding)",
#   25   │           "A1BG antisense RNA 1 (non-protein coding)"
#   26   │         ],
#   27   │         "refseq_accession": [
#   28   │           "NR_015380"
#   29   │         ],
#   30   │         "locus_type": "RNA, long non-coding",
#   31   │         "agr": "HGNC:37133",
#   32   │         "hgnc_id": "HGNC:37133",
#   33   │         "ensembl_gene_id": "ENSG00000268895",
#   34   │         "entrez_id": "503538",
#   35   │         "gene_group": [
#   36   │           "Antisense RNAs"
#   37   │         ],
#   38   │         "symbol": "A1BG-AS1",
#   39   │         "date_name_changed": "2012-08-15",
#   40   │         "location": "19q13.43",
#   41   │         "lncipedia": "A1BG-AS1",
#   42   │         "name": "A1BG antisense RNA 1",
#   43   │         "date_modified": "2013-06-27",
#   44   │         "ucsc_id": "uc002qse.3",
#   45   │         "prev_symbol": [
#   46   │           "NCRNA00181",
#   47   │           "A1BGAS",
#   48   │           "A1BG-AS"
#   49   │         ],
#   50   │         "ena": [
#   51   │           "BC040926"
#   52   │         ],
#   53   │         "gene_group_id": [
#   54   │           1987
#   55   │         ],
#   56   │         "date_symbol_changed": "2010-11-25",
#   57   │         "location_sortable": "19q13.43"
#   58   │       },

@attr.s()
class RawEntry:
    agr_id = attr.ib(validator=is_a(str))
    symbol = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str)) 
    hgnc_id = attr.ib(validator=is_a(str))
    ucsc_id = attr.ib(validator=is_a(str))
    lncipedia_id = attr.ib(validator=optional(is_a(str)))
    rnacentral_id = attr.ib(validator=optional(is_a(str)))
    previous_names: ty.List[str] = attr.ib(validator=is_a(list))
    refseq_ids: ty.List[str] = attr.ib(validator=is_a(list))
    ena_ids: ty.List[str] = attr.ib(validator=is_a(list))

    @property
    def gtrnadb_id(self):
        if self.locus_type != 'RNA, transfer':
            return None
        return ''
