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

import typing as ty

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.rnacentral.precompute.data.accession import Accession
from rnacentral_pipeline.rnacentral.precompute.data.coordinate import Coordinate
from rnacentral_pipeline.rnacentral.precompute.data.r2dt import R2dtHit
from rnacentral_pipeline.rnacentral.precompute.data.rfam import RfamHit


def partioned_accessions(so_tree, all_accessions):
    """
    Parition the list of accessions between those athat are active and those
    that are inactive. This returns a tuple of active, inactive Accession
    objects.
    """

    accessions = []
    inactive_accessions = []
    for accession in all_accessions:
        if accession["is_active"]:
            accessions.append(Accession.build(so_tree, accession))
        else:
            inactive_accessions.append(Accession.build(so_tree, accession))
    return accessions, inactive_accessions


def fix_hgnc_data(accessions):
    hgnc_entries = []
    ref_seq_entry = []
    result = []
    for accession in accessions:
        if accession.database == "hgnc":
            hgnc_entries.append(accession)
        else:
            result.append(accession)
            if accession.database == "refseq":
                ref_seq_entry.append(accession)

    if not hgnc_entries or len(ref_seq_entry) != 1:
        return accessions

    for accession in hgnc_entries:
        result.append(
            attr.evolve(
                accession,
                feature_name=ref_seq_entry[0].feature_name,
                ncrna_class=ref_seq_entry[0].ncrna_class,
            )
        )
    return result


@attr.s(frozen=True)
class Sequence:
    """
    This is the representation of sequences for the purpose of our precompute
    step.
    """

    upi = attr.ib(validator=is_a(str))
    taxid = attr.ib(validator=is_a(int))
    length = attr.ib(validator=is_a(int))
    accessions: ty.List[Accession] = attr.ib(validator=is_a(list))
    inactive_accessions: ty.List[Accession] = attr.ib(validator=is_a(list))
    is_active = attr.ib(validator=is_a(bool))
    previous_update: ty.Dict[str, str] = attr.ib(validator=is_a(dict))
    rfam_hits: ty.List[RfamHit] = attr.ib(validator=is_a(list))
    coordinates: ty.List[Coordinate] = attr.ib(validator=is_a(list))
    last_release = attr.ib(validator=is_a(int))
    r2dt_hits: ty.List[R2dtHit] = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, so_tree, data) -> "Sequence":
        """
        Given a dictonary of result from the precompute query this will build a
        SpeciesSequence.
        """

        active, inactive = partioned_accessions(so_tree, data["accessions"])
        active = fix_hgnc_data(active)
        coords = [Coordinate.build(c) for c in data["coordinates"] if c["assembly_id"]]
        previous = {}
        if data["previous"] is not None:
            previous = data["previous"]

        return cls(
            upi=data["upi"],
            taxid=data["taxid"],
            length=data["length"],
            accessions=active,
            inactive_accessions=inactive,
            is_active=not data["deleted"],
            previous_update=previous,
            rfam_hits=list(
                {
                    RfamHit.build(so_tree, r)
                    for r in data.get("rfam_hits", [])
                    if r["rfam_hit_id"]
                }
            ),
            coordinates=coords,
            last_release=data["last_release"],
            r2dt_hits=list(
                {
                    R2dtHit.build(so_tree, r)
                    for r in data.get("r2dt_hits", [])
                    if r["model_id"]
                }
            ),
        )

    @property
    def rna_id(self) -> str:
        """
        Build the RNA id which is {upi}_{taxid}.
        """

        if self.taxid is not None:
            return f"{self.upi}_{self.taxid}"
        return self.upi

    @property
    def has_coordinates(self) -> bool:
        return bool(self.coordinates)

    @property
    def chromosomes(self) -> ty.Set[str]:
        return {coord.chromosome for coord in self.coordinates}

    def is_mitochondrial(self) -> bool:
        """
        Check if this accession is mitochrondrial.
        """
        return any(a.is_mitochondrial() for a in self.accessions) or (
            self.taxid == 9606 and "MT" in self.chromosomes
        )

    def is_chloroplast(self) -> bool:
        return (
            any(a.is_chloroplast() for a in self.accessions)
            or self.taxid == 9606
            and "MT" in self.chromosomes
        )

    def species(self) -> ty.Set[str]:
        all_species: ty.Set[str] = set()
        for accession in self.accessions:
            if not accession.species:
                continue
            for species in accession.species:
                if species:
                    all_species.add(species)
        return all_species

    def domains(self) -> ty.Set[str]:
        """
        Get the set of all domains assigned to this sequence. This will ignore
        any invalid assignments (environmental samples, uncultured, or
        synthetic).
        """
        return {a.domain for a in self.accessions if a.domain}

    def has_unique_hit(self) -> bool:
        """
        Check if there is only one Rfam hit to this sequence.
        """
        return len(self.rfam_hits) == 1

    def has_r2dt_match(self) -> bool:
        return bool(self.r2dt_hits)

    def has_rfam_hit(self) -> bool:
        return bool(self.rfam_hits)
