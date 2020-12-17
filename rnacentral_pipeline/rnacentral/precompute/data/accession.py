# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

from collections import Counter

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


@attr.s(frozen=True)
class Accession:
    """
    This represents the minimal amount of information we need to produce a good
    name from the accession level data for a sequence.
    """

    gene = attr.ib(validator=optional(is_a(str)), converter=str)
    optional_id = attr.ib(validator=optional(is_a(str)))
    pretty_database = attr.ib(validator=is_a(str), converter=str)
    feature_name = attr.ib(validator=is_a(str), converter=str)
    ncrna_class = attr.ib(validator=optional(is_a(str)))
    species = attr.ib(validator=optional(is_a(str)))
    common_name = attr.ib(validator=optional(is_a(str)))
    description = attr.ib(validator=is_a(str), converter=str)
    locus_tag = attr.ib(validator=optional(is_a(str)))
    organelle = attr.ib(validator=optional(is_a(str)))
    lineage = attr.ib(validator=optional(is_a(str)))
    all_species = attr.ib(validator=is_a(tuple), converter=tuple)
    all_common_names = attr.ib(validator=is_a(tuple), converter=tuple)
    so_rna_type = attr.ib(validator=is_a(str), converter=str)

    @classmethod
    def build(cls, data):
        """
        Create a new Accession from the given dict. This assumes the dict has
        keys with the same names as accession fields.
        """
        return cls(**data)

    @property
    def database(self):
        """
        The normalized (lowercase) database name.
        """
        return self.pretty_database.lower()

    @property
    def rna_type(self):
        """
        Get a single INSDC RNA type for this accession.
        """

        if self.feature_name == "ncRNA":
            return self.ncrna_class
        return self.feature_name

    @property
    def domain(self):
        """
        Get the domain, if any, that is assigned to this accession. This will
        not include uncultured, environmental or synthetic domains.
        """

        if self.lineage.startswith("other sequences"):
            return None

        parts = [p.strip() for p in self.lineage.split(";")]
        domain = parts[0]
        if domain == "cellular organisms":
            domain = parts[1]

        if (
            "uncultured" in domain
            or "environmental" in domain
            or "synthetic" in domain
            or "artificial" in domain
            or "unclassified" in domain
        ):
            return None

        return domain

    @property
    def masked_description(self):
        """
        Compute a masked description. This will do things like strip out
        '10-mer' and such. The description returned is suitable for entropy
        computation, but as the description that is displayed to the user.
        """

        raw = self.description.lower()  # pylint: disable=no-member
        allowed = set(string.ascii_lowercase + string.digits + " ")
        counts = Counter(r for r in raw if r in allowed)
        rep = counts.most_common(1)[0][0]
        masked = re.sub(r"(\d+-mer)", lambda m: rep * len(m.groups(0)[0]), raw)
        masked = re.sub(r"5'-(.+)-3'", "", masked)
        masked = "".join(m for m in masked if m in allowed)
        masked = re.sub(r"\s+", " ", masked)
        return masked

    def is_mitochondrial(self):
        """
        Check if this accession is mitochrondrial.
        """
        return "mitochondri" in self.description or (
            self.organelle and "mitochondri" in self.organelle
        )

    def is_chloroplast(self):
        return "chloroplast" in self.description or (
            self.organelle and "chloroplast" in self.organelle
        )
