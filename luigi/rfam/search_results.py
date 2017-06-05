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

import pandas as pd

INFERNAL_HEADER = {0: "hit_num",
                   1: "rfam_id",
                   2: "rfam_acc",
                   3: "urs",
                   4: "accession",
                   5: "clan_acc",
                   6: "mdl",
                   7: "mdl_from",
                   8: "mdl_to",
                   9: "seq_from",
                   10: "seq_to",
                   11: "strand",
                   12: "trunc",
                   13: "pass",
                   14: "gc",
                   15: "bias",
                   16: "score",
                   17: "e_value",
                   18: "inc",
                   19: "olp"}


def parse(filename):
    # read with pandas df, variable space sep, only first 20 columns
    df_tbl = pd.read_table(filename,
                           comment="#",
                           engine="python",
                           sep=r"\s*|-\n",
                           header=None,
                           usecols=range(20))

    df_tbl.columns = [INFERNAL_HEADER[column] for column in df_tbl]
    rel_columns = [3, 0, 2, 5, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19]
    return df_tbl[rel_columns]
