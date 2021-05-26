# -*- coding: utf-8 -*-

from __future__ import annotations

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import enum
import re
import typing as ty
from pathlib import Path

from Bio import SeqIO

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline.databases.data import RibovoreResult

TRNAS = {
    "SO:0000254",
    "SO:0000259",
    "SO:0000268",
    "SO:0000260",
    "SO:0000266",
    "SO:0000256",
    "SO:0000270",
    "SO:0000261",
    "SO:0000273",
    "SO:0000272",
    "SO:0000258",
    "SO:0000263",
    "SO:0000269",
    "SO:0000264",
    "SO:0000271",
    "SO:0005857",
    "SO:0000766",
    "SO:0000265",
    "SO:0000257",
    "SO:0001036",
    "SO:0000262",
    "SO:0001172",
    "SO:0002129",
    "SO:0000267",
}


@enum.unique
class Source(enum.Enum):
    crw = enum.auto()
    ribovision = enum.auto()
    rfam = enum.auto()
    rnase_p = enum.auto()
    gtrnadb = enum.auto()

    @classmethod
    def build(cls, name: ty.Union[str, Source]) -> Source:
        if isinstance(name, Source):
            return name
        name = name.lower()
        if hasattr(cls, name):
            return getattr(cls, name)
        if name == 'rnase p database':
            return Source.rnase_p
        raise ValueError(f"Unknown database name {name}")

    def result_directory(self) -> str:
        if self is Source.crw:
            return "crw"
        if self is Source.ribovision:
            return "ribovision"
        if self is Source.rfam:
            return "rfam"
        if self is Source.rnase_p:
            return "rnasep"
        if self is Source.gtrnadb:
            return "gtrnadb"
        raise ValueError(f"Could not find results for {self}")


@attr.s()
class ModelInfo(object):
    model_name: str = attr.ib(validator=is_a(str))
    is_intronic: bool = attr.ib(validator=is_a(bool))
    so_term: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int))
    accessions: ty.List[str] = attr.ib(validator=is_a(list))
    source: Source = attr.ib(validator=is_a(Source))
    length: ty.Optional[int] = attr.ib(validator=optional(is_a(int)))
    basepairs = attr.ib(validator=optional(is_a(int)))
    cell_location: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))

    @property
    def rna_type(self):
        if self.so_term in {"SO:0000650", "SO:0000651", "SO:0000652", "SO:0001001"}:
            return "rRNA"
        if self.so_term in {"SO:0000587", "SO:0000603"}:
            return "autocatalytically_spliced_intron"
        if self.so_term in TRNAS:
            return "tRNA"
        if self.so_term == "SO:0000386":
            return "RNase_P_RNA"
        raise ValueError("No RNA type for: %s" % self)

    def writeable(self):
        return [
            self.model_name,
            self.taxid,
            self.rna_type,
            self.so_term,
            self.cell_location,
            self.source.name,
            self.length,
        ]

    def fraction_paired(self) -> ty.Optional[float]:
        if self.length is not None and self.basepairs is not None:
            return float(self.basepairs) / float(self.length)
        return None


@attr.s(hash=True)
class ModelDatabaseInfo:
    name = attr.ib(validator=is_a(str))
    db_id = attr.ib(validator=is_a(int))
    source = attr.ib(validator=is_a(Source))
    alias = attr.ib(validator=optional(is_a(str)))
    length = attr.ib(validator=optional(is_a(int)))
    basepairs = attr.ib(validator=optional(is_a(int)))

    @classmethod
    def build(cls, raw) -> ModelDatabaseInfo:
        return cls(
            name=raw["model_name"],
            db_id=raw["model_id"],
            source=getattr(Source, raw["model_source"]),
            alias=raw["model_alias"],
            length=raw["model_length"],
            basepairs=raw["model_basepairs"],
        )


@attr.s(hash=True)
class R2DTResultInfo(object):
    urs = attr.ib(validator=is_a(str))
    db_info = attr.ib(validator=is_a(ModelDatabaseInfo))
    source = attr.ib(validator=is_a(Source))
    path = attr.ib(validator=is_a(Path))

    @property
    def model_name(self):
        return self.db_info.name

    @property
    def model_db_id(self):
        return self.db_info.db_id

    @property
    def model_alias(self):
        return self.db_info.alias

    @property
    def model_basepairs(self) -> int:
        if self.db_info.basepairs is None:
            raise ValueError(f"No model length for {self.db_info}")
        return self.db_info.basepairs

    @property
    def model_length(self) -> int:
        if self.db_info.length is None:
            raise ValueError(f"No model length for {self.db_info}")
        return self.db_info.length

    @property
    def svg(self) -> Path:
        base = self.path / "svg"
        paths = list(base.glob(f"{self.urs}*.svg"))
        if not paths:
            raise ValueError(f"Could not figure out svg filename for {self}")
        if len(paths) > 1:
            raise ValueError("Too many possible svg files")
        return paths[0]

    @property
    def fasta(self) -> Path:
        return self.path / "fasta" / self.__filename__("fasta")

    @property
    def source_directory(self) -> Path:
        base = (self.path / "..").resolve()
        if self.source == Source.ribovision:
            parts = self.model_name.split("_", 3)
            if parts[1] == "LSU":
                return base / "ribovision-lsu"
            elif parts[1] == "SSU":
                return base / "ribovision-ssu"
            raise ValueError("Could not find correct data path: %s" % parts)

        if self.source == Source.rfam and self.model_name in {"RF00005", "tRNA"}:
            return base / "RF00005"
        return base / self.source.result_directory()

    @property
    def overlaps(self) -> Path:
        base = self.source_directory
        paths = list(base.glob(f"{self.urs}*.overlaps"))
        if not paths:
            print(self)
            raise ValueError(
                f"Could not figure out overlaps filename for {self.source_directory}/{self.urs}")
        if len(paths) > 1:
            raise ValueError(f"Too many overlaps files for {self.urs}")
        return paths[0]

    def publish_path(self, suffix="", compressed=False) -> Path:
        publish = Path(self.urs[0:3])
        for start in range(4, 11, 2):
            publish = publish / self.urs[start: (start + 2)]
        append = ""
        if suffix:
            append = f"-{suffix}"
        extension = "svg"
        if compressed:
            extension += ".gz"
        return publish / f"{self.urs}{append}.{extension}"

    def validate(self):
        assert self.svg.exists(), "Missing SVG file (%s) for %s" % (self.svg, self)
        assert self.fasta.exists(), "Missing FASTA file (%s) for %s" % (
            self.fasta,
            self,
        )
        assert self.overlaps.exists(), "Missing overlaps (%s) for %s" % (
            self.overlaps,
            self,
        )
        assert self.source_directory.exists(), "Missing source (%s) for %s" % (
            self.source_directory,
            self,
        )

    def has_ribovore(self):
        if self.source in {Source.crw, Source.ribovision}:
            return True
        if self.source == Source.rfam and self.model_name not in {"RF00005", "tRNA"}:
            return True
        return False

    def has_hit_info(self):
        return self.has_ribovore()

    def __filename__(self, extension):
        if self.source == Source.gtrnadb and extension == "fasta":
            return f"{self.urs}.{extension}"
        if self.source == Source.rfam and not self.model_name.startswith("RF"):
            if extension == "fasta":
                return f"{self.urs}.{extension}"
            assert self.model_alias.startswith(
                "RF"), f"No existing alias for {self}"
            return f"{self.urs}-{self.model_alias}.{extension}"
        return f"{self.urs}-{self.model_name}.{extension}"


@attr.s()
class R2DTResult(object):
    info = attr.ib(validator=is_a(R2DTResultInfo))
    hit_info = attr.ib(validator=optional(is_a(RibovoreResult)), default=None)

    @classmethod
    def from_info(cls, info: R2DTResultInfo, hit_info=None):
        return cls(info, hit_info=hit_info)

    @property
    def urs(self):
        return self.info.urs

    @property
    def model_id(self):
        return self.info.model_db_id

    @property
    def source(self):
        return self.info.source

    @property
    def publish_path(self):
        return self.info.publish_path

    def svg(self):
        """
        Process a single SVG file into the requried data. This produce an array
        that can be written to CSV for import into the database.
        """
        with self.info.svg.open("r") as raw:
            return raw.read().replace("\n", "")

    def dot_bracket(self):
        """
        Parse the fasta/dot_bracket file for the given urs and model to extract the
        dot_bracket string for the pairing of the URS in that model. This assumes
        that there is only a single record in the file and it has a sequence then
        dot_bracket string which are the same length.
        """
        with self.info.fasta.open("r") as raw:
            record = SeqIO.read(raw, "fasta")
            seq_dot = str(record.seq)
            sequence = re.match(r"^(\w+)", seq_dot).group(1)
            dot_bracket = re.sub(r"^\w+", "", seq_dot)
            assert len(sequence) == len(dot_bracket)
            return dot_bracket

    def basepair_count(self):
        return self.dot_bracket().count("(")

    def modeled_length(self):
        return len(self.dot_bracket())

    def overlap_count(self):
        with self.info.overlaps.open("r") as raw:
            return int(raw.readline().strip())

    @property
    def model_basepairs(self):
        return self.info.model_basepairs

    @property
    def model_length(self):
        return self.info.model_length

    def writeable(self):
        model_start = None if not self.hit_info else self.hit_info.mfrom
        model_stop = None if not self.hit_info else self.hit_info.mto
        sequence_start = None if not self.hit_info else self.hit_info.bfrom
        sequence_stop = None if not self.hit_info else self.hit_info.bto
        sequence_coverage = None if not self.hit_info else self.hit_info.bcov

        return [
            self.urs,
            self.model_id,
            self.dot_bracket(),
            self.overlap_count(),
            self.basepair_count(),
            model_start,
            model_stop,
            sequence_start,
            sequence_stop,
            sequence_coverage,
            True,
        ]


@attr.s()
class ShowInfo:
    urs = attr.ib(validator=is_a(str))
    model_id = attr.ib(validator=is_a(int))
    model_length = attr.ib(validator=is_a(int))
    model_basepairs = attr.ib(validator=is_a(int))
    sequence_length = attr.ib(validator=is_a(int))
    modeled_length = attr.ib(validator=is_a(int))
    modeled_basepairs = attr.ib(validator=is_a(int))

    @classmethod
    def from_raw(cls, raw) -> ShowInfo:
        return cls(
            urs=raw["urs"],
            model_id=raw["model_id"],
            model_length=raw["model_length"],
            model_basepairs=raw["model_basepairs"],
            sequence_length=raw["sequence_length"],
            modeled_length=raw["modeled_length"],
            modeled_basepairs=raw["modeled_basepairs"],
        )

    @classmethod
    def from_result(cls, result: R2DTResult) -> ShowInfo:
        assert result.model_length is not None, "Cannot show without length"
        assert result.model_basepairs is not None, "Cannot show without bps"

        return cls(
            urs=result.urs,
            model_id=result.model_id,
            model_length=result.model_length,
            model_basepairs=result.model_basepairs,
            sequence_length=result.sequence_length,
            modeled_length=result.modeled_length(),
            modeled_basepairs=result.basepair_count(),
        )

    def observed_paired(self) -> float:
        return float(self.modeled_basepairs) / self.modeled_length

    def model_paired(self) -> float:
        return float(self.model_basepairs) / self.model_length

    def showable(self) -> bool:
        observed = round(self.observed_paired(), 1)
        modeled = round(self.model_paired(), 1)
        return observed >= modeled

    def writeable(self) -> ty.List[str]:
        return [self.urs, str(self.model_id), str(self.showable())]
