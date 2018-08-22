# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import operator as op

from rnacentral_pipeline.writers import MultiCsvOutput


def write_annotations(parser, *args, **kwargs):
    writer = MultiCsvOutput.build(
        parser,
        annotations={
            'transformer': op.methodcaller('writeable'),
        },
        publications={
            'transformer': op.methodcaller('writeable_publications'),
        },
        publication_mappings={
            'transformer': op.methodcaller('writeable_publication_mappings'),
        },
        terms={
            'transformer': op.methodcaller('writeable_ontology_terms'),
        },
    )
    writer(*args, **kwargs)
