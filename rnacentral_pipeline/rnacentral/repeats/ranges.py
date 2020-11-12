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

import csv
import pickle
import subprocess as sp
import typing as ty
from pathlib import Path

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.ensembl.metadata import databases as db


@attr.s(slots=True)
class Info:
    assembly_id = attr.ib(validator=is_a(str))
    compressed = attr.ib(validator=is_a(Path), converter=Path)
    index_file = attr.ib(validator=is_a(Path), converter=Path)
    chromosome_column = attr.ib(validator=is_a(int))
    start_column = attr.ib(validator=is_a(int))
    stop_column = attr.ib(validator=is_a(int))

    @classmethod
    def from_file(cls, path: Path) -> "Info":
        if not path.exists():
            raise ValueError(f"Expected info path {path} does not exist")
        with path.open("rb") as raw:
            info = pickle.load(raw)
            info.validate()
            return info

    @classmethod
    def from_directory(cls, path: Path) -> "Info":
        return cls.from_file(path / "info.pickle")

    def validate(self):
        assert self.start_column != self.stop_column != self.chromosome_column
        assert (
            self.compressed.exists()
        ), f"Compressed file {self.compressed} does not exist"
        assert (
            self.index_file.exists()
        ), f"Compressed file {self.index_file} does not exist"

    def dump(self):
        self.validate()
        with open("info.pickle", "wb") as out:
            pickle.dump(data, out)


@attr.s(frozen=True, slots=True)
class RepeatRanges:
    """
    This represents all repeat ranges in an assembly
    """

    info = attr.ib(validator=is_a(Info))

    @classmethod
    def from_info(cls, info: Info) -> "RepeatRanges":
        """
        Load the RepeatRanges from the given path.
        """
        info.validate()
        return cls(info=info)

    def overlaps(self, chromosome: str, start: int, stop: int) -> ty.Iterable[ty.List]:
        """
        Fetch all intervals that overlap the given start/stop in the given
        chromosome.
        """
        query = "{chromosome}:{start}-{stop}"
        process = sp.Popen(["tabix", str(self.info.compressed), query])
        for line in process.stdout:
            yield line.strip().split()

    def enveloped_by(
        self, chromosome: str, start: int, stop: int
    ) -> ty.Iterable[ty.List]:
        """
        Find all intervals which envelop the given start/stop.
        """

        for overlap in self.overlaps(chromosome, start, stop):
            found_start = overlap[self.info.start_column]
            found_stop = overlap[self.info.stop_column]
            if found_start <= start and found_stop >= found_stop:
                yield overlap

    def is_enveloped(self, chromosome: str, start: int, stop: int) -> bool:
        intervals = self.enveloped_by(chromosome, start, stop)
        return bool(next(intervals, False))


def build_bed_directory(
    assembly: str, path: Path, chromosome_column=1, start_column=2, stop_column=3
):
    assert path.is_dir(), "{path} must be an existing directory"
    compressed = path / f"{assembly}.bed.bgz"
    indexed = path / f"{assembly}.bed.bgz.tbi"
    assert compressed.is_file(), f"Compressed file {compressed} must exist"
    assert indexed.is_file(), f"Indexed file {indexed} must exist"
    info = Info(
        assembly_id=assembly,
        compressed=compressed,
        index_file=indexed,
        chromosome_column=chromosome_column,
        start_column=start_column,
        stop_column=stop_column,
    )
    info.dump()


def find_databases(connection_handle, assembly_file):
    assemblies = {}
    for (ensembl_url, assembly_id, division) in csv.reader(assembly_file):
        assemblies[ensembl_url] = {
            "ensembl_url": ensembl_url,
            "assembly_id": assembly_id,
            "division": division,
        }
    specs = json.load(connection_handle)
    for conn_name, spec in specs.items():
        if not conn_name.lower().startswith("ensembl"):
            continue
        connection = pymysql.Connection(**spec)
        seen = set()
        for db_name in db.databases(connection):
            for key, info in assemblies.items():
                if db_name.startswith(key):
                    if info["assembly_id"] in seen:
                        raise ValueError(f"Dupcliateion assembly {info['assembly_id']}")
                    yield [info["assembly_id"], conn_name, db_name]
                    seen.add(info["assembly_id"])
                    break
            else:
                LOGGER.warn("Could not find assembly for %s", db_name)


# def from_ensembl_fasta(assembly: str, path: Path) -> RepeatRanges:
#     """
#     Parse a fasta file from Ensembl that contains soft masked repeat information
#     and produce a RepeatRanges representing that file. This assumes that the
#     id for each entry in the file is the chromosome name.
#     """
#     ranges = RepeatRanges(assembly=assembly)
#     for record in SeqIO.parse(path, "fasta"):
#         nts = enumerate(record.seq, start=1)
#         masked = more.split_at(nts, lambda nt: nt[1].upper() == nt[1])
#         filtered = filter(None, masked)
#         compressed = ((m[0][0], m[-1][0]) for m in filtered)
#         for (start, stop) in compressed:
#             ranges.insert(record.id, start, stop)
#     return ranges
