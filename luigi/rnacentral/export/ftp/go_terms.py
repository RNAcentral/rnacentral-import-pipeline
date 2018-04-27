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
import itertools as it

from rnacentral.psql import PsqlWrapper


QUERY = """
select
    pre.id,
    pre.upi,
    pre.taxid,
    go_terms.ontology_term_id,
    string_agg('Rfam:' || hits.rfam_model_id, '|') as models
from rnc_rna_precomputed pre
join rfam_model_hits hits on hits.upi = pre.upi
join rfam_go_terms go_terms on hits.rfam_model_id = go_terms.rfam_model_id
where
    exists(select 1 from qa_status qa where qa.upi = pre.upi and qa.taxid = pre.taxid and qa.has_issue = false)
    and exists(select 1 from xref where xref.upi = pre.upi and xref.taxid = pre.taxid and xref.deleted = 'N')
group by pre.id, go_terms.ontology_term_id
"""

UNCULTURED_TAXIDS = {
    155900,
    198431,
    415540,
    358574,
    81726,
}


def annotations(config):
    """
    Get the list of all annotations implied by Rfam matches.
    """

    psql = PsqlWrapper(config)
    for annotation in psql.copy_to_iterable(QUERY):
        yield annotation


def exclude_uncultured(annotation):
    """
    Detect if the annotation is from an uncultured organism and thus should be
    excluded.
    """
    return annotation['taxid'] in UNCULTURED_TAXIDS


def valid_annotations(config):
    """
    Filter all annotations to select only the valid ones.
    """

    data = annotations(config)
    return it.ifilter(exclude_uncultured, data)


def export(config, handle):
    """
    Write the Rfam based GO annotations to the given handle.
    """

    writer = csv.DictWriter(
        handle,
        ['id', 'ontology_term_id', 'models'],
        extrasaction='ignore',
    )
    writer.writerows(valid_annotations(config))
