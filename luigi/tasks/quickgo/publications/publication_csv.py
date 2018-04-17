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


class QuickGoPublicationCsv(CsvWriter):
    headers = [
        'ref_pubmed_id',
        'authors',
        'location',
        'title',
    ]

    def requires(self):
        conf = quickgo()
        return [
            FetchTask(remote_path=conf.data_file, local_path=conf.annotations),
        ]

    def data(self):
        filename = self.requires()[0].output().fn
        with gzip.open(filename, 'r') as raw:
            seen = set()
            for go_term in parser(raw):
                publication = go_term.writable_publication()
                if publication not in seen:
                    yield publication
