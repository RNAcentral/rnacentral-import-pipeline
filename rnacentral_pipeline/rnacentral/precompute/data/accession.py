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

import re
import string
import typing as ty
from collections import Counter

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline.databases.data import RnaType
from rnacentral_pipeline.databases.data import Database


@attr.s(frozen=True)
class Accession:
    """
    This represents the minimal amount of information we need to produce a good
    name from the accession level data for a sequence.
    """

    gene = attr.ib(validator=optional(is_a(str)))
    optional_id = attr.ib(validator=optional(is_a(str)))
    database = attr.ib(validator=is_a(Database))
    species = attr.ib(validator=optional(is_a(str)))
    common_name = attr.ib(validator=optional(is_a(str)))
    description = attr.ib(validator=is_a(str), converter=str)
    locus_tag = attr.ib(validator=optional(is_a(str)))
    organelle = attr.ib(validator=optional(is_a(str)))
    lineage = attr.ib(validator=optional(is_a(str)))
    all_species = attr.ib(validator=is_a(tuple), converter=tuple)
    all_common_names = attr.ib(validator=is_a(tuple), converter=tuple)
    rna_type = attr.ib(validator=is_a(RnaType))
    is_active = attr.ib(validator=is_a(bool))

    @classmethod
    def build(cls, so_tree, data) -> "Accession":
        """
        Create a new Accession from the given dict. This assumes the dict has
        keys with the same names as accession fields.
        """
        rna_type = None
        if data["so_rna_type"]:
            rna_type = RnaType.from_so_term(so_tree, data["so_rna_type"])
        else:
            name = data["ncrna_class"] or data["feature_name"]
            rna_type = RnaType.from_insdc_term(so_tree, name)
        return cls(
            gene=data["gene"],
            optional_id=data["optional_id"],
            database=Database.build(data["database"]),
            species=data["species"],
            common_name=data["common_name"],
            description=data["description"],
            locus_tag=data["locus_tag"],
            organelle=data["organelle"],
            lineage=data["lineage"],
            all_species=tuple(data["all_species"]),
            all_common_names=tuple(data["all_common_names"]),
            rna_type=rna_type,
            is_active=data["is_active"],
        )

    @property
    def domain(self) -> ty.Optional[str]:
        """
        Get the domain, if any, that is assigned to this accession. This will
        not include uncultured, environmental or synthetic domains.
        """

        if not self.lineage:
            return None
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
            or "other entries" in domain
        ):
            return None

        return domain

    @property
    def masked_description(self) -> str:
        """
        Compute a masked description. This will do things like strip out
        '10-mer' and such. The description returned is suitable for entropy
        computation, but as the description that is displayed to the user.
        """

        raw = self.description.lower()
        allowed = set(string.ascii_lowercase + string.digits + " ")
        counts = Counter(r for r in raw if r in allowed)
        rep = counts.most_common(1)[0][0]
        masked = re.sub(r"(\d+-mer)", lambda m: rep * len(m.group(1)[0]), raw)
        masked = re.sub(r"5'-(.+)-3'", "", masked)
        masked = "".join(m for m in masked if m in allowed)
        masked = re.sub(r"\s+", " ", masked)
        return masked

    @property
    def pretty_database(self):
        return self.database.pretty()

    def is_mitochondrial(self) -> bool:
        """
        Check if this accession is mitochrondrial.
        """
        found = "mitochondri" in self.description or (
            self.organelle and "mitochondri" in self.organelle
        )
        return bool(found)

    def is_chloroplast(self) -> bool:
        found = "chloroplast" in self.description or (
            self.organelle and "chloroplast" in self.organelle
        )
        return bool(found)
