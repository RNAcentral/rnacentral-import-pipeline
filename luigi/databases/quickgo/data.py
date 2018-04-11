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

import attr
import requests
from functools32 import lru_cache


ANN_URL = 'http://www.ebi.ac.uk/QuickGO/annotations?geneProductId={upi}'
TERM_URL = 'http://www.ebi.ac.uk/QuickGO/term/{go_id}'
INFO_URL = 'http://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}'


@attr.s()
class GoTerm(object):
    """
    This represents what
    """

    go_id = attr.ib()
    name = attr.ib()
    definition = attr.ib()

    @lru_cache
    @classmethod
    def from_id(cls, go_id):
        """
        Will create a GoTerm with the associated metadata by fetching information
        using the given GO ID.
        """

        assert go_id
        response = requests.get(INFO_URL.format(go_id=go_id))
        response.raise_for_status()
        data = response.json()
        assert data["numberOfHits"] == 1, "Non-unique ID, somehow"
        data = data['results'][0]
        assert data['id'] == go_id, "Got info about wrong term, somehow"
        return cls(
            go_id=go_id,
            name=data['name'],
            definition=data['definition']['text'],
        )

    @property
    def url(self):
        """
        Get the url for the human readable page about this GO term.
        """
        return TERM_URL.format(go_id=self.go_id)


@attr.s()
class GoTermAnnotation(object):
    """
    This is the simpliest possible representation of AmiGO annotations about
    RNAcentral sequences.
    """

    upi = attr.ib(validator=is_a(basestring))
    qualifier = attr.ib(validator=is_a(basestring))
    go_id = attr.ib(validator=is_a(basestring))
    publications = attr.ib(validator=is_a(list))

    @property
    def url(self):
        """
        The URL for this GO annotation.
        """
        return ANN_URL.format(upi=self.upi)

    def go_term(self):
        """
        Fetch the information about the GoTerm this annotation uses.
        """
        return GoTerm.from_id(self.go_id)
