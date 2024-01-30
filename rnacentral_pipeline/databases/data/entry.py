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

import itertools as it
import json
import logging
import operator as op
import typing as ty
from collections import Counter

import attr
from attr.validators import and_
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.helpers.hashes import crc64, md5

from . import utils
from .features import SequenceFeature
from .go_annotations import GoTermAnnotation
from .references import IdReference, Reference
from .regions import Exon, SequenceRegion
from .secondary_structure import SecondaryStructure

LOGGER = logging.getLogger(__name__)

FEATURE_TYPE_RNAS = set(
    [
        "SO:0000252",
        "SO:0000253",
        # 'precursor_RNA',
        "SO:0000584",
        "SO:0000673",
    ]
)


@attr.s(frozen=True)
class Entry:
    """
    This represents an RNAcentral entry that will be imported into the
    database. It should contain all the information needed to define all the
    data that is loaded from expert databases. For example it should contain
    information for rna (sequence), rnc_accessions, rnc_coordinates, and so
    forth.
    """

    # Also known as external_id
    primary_id = attr.ib(validator=is_a(str), converter=str)
    accession = attr.ib(validator=is_a(str), converter=str)
    ncbi_tax_id = attr.ib(validator=is_a(int))
    database = attr.ib(
        validator=is_a(str),
        converter=lambda s: str(s.upper()),
    )
    sequence = attr.ib(validator=is_a(str), converter=str)
    regions: ty.List[SequenceRegion] = attr.ib(validator=is_a(list))
    rna_type = attr.ib(
        validator=utils.matches_pattern(utils.SO_PATTERN),
        converter=utils.as_so_term,
    )
    url = attr.ib(validator=is_a(str), converter=str)
    seq_version = attr.ib(
        validator=and_(is_a(str), utils.matches_pattern(r"^\d+$")),
        converter=lambda r: str(int(float(r))),
    )

    note_data = utils.possibly_empty(dict)
    xref_data = utils.possibly_empty(dict)

    related_sequences = utils.possibly_empty(list)
    related_diseases = utils.possibly_empty(list)

    chromosome: str = utils.optionally(str)
    species: str = utils.optionally(str)
    common_name: str = utils.optionally(str)
    lineage: str = utils.optionally(str)
    gene: str = utils.optionally(str)
    locus_tag: str = utils.optionally(str)
    optional_id: str = utils.optionally(str)
    product: str = utils.optionally(str)
    parent_accession: str = utils.optionally(str)
    non_coding_id: str = utils.optionally(str)
    project: str = utils.optionally(str)
    keywords: str = utils.optionally(str)
    organelle: str = utils.optionally(str)
    anticodon: str = utils.optionally(str)
    experiment: str = utils.optionally(str)
    function: str = utils.optionally(str)
    inference: str = utils.optionally(str)
    standard_name: str = utils.optionally(str)
    description: str = utils.optionally(str)
    mol_type: str = utils.optionally(str)
    is_composite: str = utils.optionally(str)

    location_start = utils.optionally(int)
    location_end = utils.optionally(int)

    gene_synonyms = utils.possibly_empty(list)
    references = utils.possibly_empty(list)

    secondary_structure = utils.possibly_empty(SecondaryStructure)
    features = utils.possibly_empty(list)
    interactions = utils.possibly_empty(list)
    go_annotations: ty.List[GoTermAnnotation] = utils.possibly_empty(list)

    @property
    def database_name(self):
        """
        Get the database name in upper case.
        """
        return self.database.upper()  # pylint: disable=E1101

    @property
    def exons(self) -> ty.List[Exon]:
        exons = []
        for region in self.regions:
            exons.extend(region.exons)
        return exons

    @property
    def db_xrefs(self) -> str:
        """
        Return a JSON encoded dict representing the xref data.
        """
        return json.dumps(self.xref_data)

    @property
    def note(self) -> str:
        """
        Return a JSON encoded dictionary representing the note data.
        """
        data = self.note_data
        if "url" not in data:
            data["url"] = self.url
        return json.dumps(data)

    @property
    def feature_name(self) -> str:
        """
        Return the feature for the RNA type.
        """
        if self.rna_type in FEATURE_TYPE_RNAS:
            return utils.SO_INSDC_MAPPING[self.rna_type]
        return "ncRNA"

    @property
    def ncrna_class(self) -> str:
        """
        The ncRNA class. If the feature type is not ncRNA this this will be the
        empty string.
        """
        if self.feature_name != "ncRNA":
            return ""
        return utils.SO_INSDC_MAPPING.get(self.rna_type, "other")

    @property
    def gene_synonym(self) -> str:
        """
        Returns a comma separated list of gene synonyms.
        """
        if self.gene_synonyms:
            return ",".join(self.gene_synonyms)
        else:
            return ""

    @property
    def feature_location_start(self):
        """
        This will compute the feature location start if it is not set,
        otherwise this will use the set one.
        """

        if self.location_start is not None:
            return self.location_start
        if not self.exons:
            return 1
        return min(e.start for e in self.exons)

    @property
    def feature_location_end(self):
        """
        This will compute the feature location end if it is not set,
        otherwise this will use the set one.
        """

        if self.location_end is not None:
            return self.location_end
        if not self.exons:
            return len(self.sequence) + 1
        return max(e.stop for e in self.exons)

    def crc64(self) -> str:
        """
        Compute a CRC64 check sum for the sequence.
        """
        return str(crc64(self.sequence))

    def md5(self) -> str:
        """
        Compute an MD5 hash of the sequence.
        """
        return str(md5(self.sequence.encode("utf-8")))

    def is_valid(self) -> bool:
        """
        Detect if this entry is valid. This means it is neither too short (< 10
        nt) not too long (> 1000000 nts) and has less than 10% N's.
        """

        assert self.description, "All entries must have a description"
        length = len(self.sequence)
        if length < 10:
            LOGGER.warn("%s is too short (%s)", self.accession, length)
            return False

        if length > 1000000:
            LOGGER.warn("%s is too long (%s)", self.accession, length)
            return False

        counts = Counter(self.sequence)
        fraction = float(counts.get("N", 0)) / float(len(self.sequence))
        if fraction > 0.1:
            LOGGER.warn(
                "%s has too many (%i/%i) N's", self.accession, counts["N"], length
            )
            return False

        if self.rna_type == "SO:0000234":
            LOGGER.warn("Skipping a mRNA")
            return False

        if "U" in self.sequence:
            raise ValueError("Sequence %s contains a U, this is very bad" % self)

        return True

    def human_rna_type(self) -> str:
        return utils.SO_INSDC_MAPPING[self.rna_type].replace("_", " ")

    def write_ac_info(self) -> ty.Iterable[ty.List[ty.Optional[str]]]:
        if not self.is_valid():
            return
        yield [
            self.accession,
            self.parent_accession,
            self.seq_version,
            self.feature_location_start,
            self.feature_location_end,
            self.feature_name,
            self.is_composite,
            self.non_coding_id,
            self.database_name,
            self.primary_id,
            self.optional_id,
            self.project,
            self.keywords,
            self.description,
            self.species,
            self.organelle,
            self.chromosome,
            self.experiment,
            self.function,
            self.gene,
            self.gene_synonym,
            self.inference,
            self.locus_tag,
            self.mol_type,
            self.ncrna_class,
            self.note,
            self.product,
            self.standard_name,
            self.db_xrefs,
            self.rna_type,
        ]

    def write_secondary_structure(self) -> ty.List[str]:
        if not self.is_valid():
            return []
        # pylint: disable=no-member
        return self.secondary_structure.writeable(self.accession)

    def write_sequence(self) -> ty.Iterable[ty.List[str]]:
        if not self.is_valid():
            return
        yield [
            self.crc64(),
            len(self.sequence),
            self.sequence,
            self.database_name,
            self.accession,
            self.optional_id,
            self.seq_version,
            self.ncbi_tax_id,
            self.md5(),
        ]

    def write_seq_short(self) -> ty.Iterable[ty.List[str]]:
        if len(self.sequence) <= 4000:
            return self.write_sequence()
        return []

    def write_seq_long(self) -> ty.Iterable[ty.List[str]]:
        if len(self.sequence) > 4000:
            return self.write_sequence()
        return []

    def write_refs(self):
        refs = filter(lambda r: isinstance(r, Reference), self.references)
        return self.__write_part__(refs)

    def write_ref_ids(self):
        refs = self.references
        refs = filter(lambda r: isinstance(r, IdReference), refs)
        return self.__write_part__(refs)

    def write_related_sequences(self):
        return self.__write_part__(self.related_sequences)

    def write_sequence_features(self):
        for feature in self.features:
            yield feature.writeable(self.accession, self.ncbi_tax_id)

        for related in self.related_sequences:
            features = related.write_features(self.accession, self.ncbi_tax_id)
            for feature in features:
                yield feature

    def write_sequence_regions(self):
        return self.__write_part__(self.regions)

    def write_interactions(self):
        for interaction in self.interactions:
            yield interaction.writeable()

    def write_ontology_terms(self) -> ty.Iterable[ty.List[str]]:
        yield [self.rna_type]

    def __write_part__(self, attribute, method_name="writeable"):
        if not self.is_valid():
            return []
        method = op.methodcaller(method_name, self.accession)
        writeable = map(method, attribute)
        return it.chain.from_iterable(writeable)
