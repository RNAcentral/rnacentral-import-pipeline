import csv
import json
import re
import typing as ty
from pathlib import Path
from urllib.parse import urlencode

import attr
import numpy as np
import pandas as pd
from attr.validators import instance_of as is_a

from rnacentral_pipeline import schemas
from rnacentral_pipeline.output_format import is_parquet
from rnacentral_pipeline.parquet_writers import parquet_writer

REDI_BASE_URL = "https://rediportal.cloud.ba.infn.it/cgi/atlas/getpage_dev.py"

REPEAT_COLUMNS = ("Repeats", "Repeat", "repeat", "Repeat_type", "type", "Location")
STRAND_COLUMNS = ("Strand", "strand")


def species_prefix(genome_build: str) -> str:
    match = re.match(r"[A-Za-z]+", genome_build)
    if not match:
        raise ValueError(f"Cannot determine species prefix for {genome_build}")
    return match.group(0)


def region_id(chromosome: str, position: ty.Any, genome_build: str) -> str:
    return f"{chromosome}_{position}_{genome_build}"


def editing_url(chromosome: str, start: int, stop: int, genome_build: str) -> str:
    query = {
        "query1": f"chr{chromosome}:{start}-{stop}",
        "query10": genome_build,
        "query9": species_prefix(genome_build),
    }
    return f"{REDI_BASE_URL}?{urlencode(query)}"


def repeat_type(raw: dict[str, ty.Any]) -> str:
    for column in REPEAT_COLUMNS:
        value = raw.get(column)
        if value and str(value).strip():
            return str(value).strip()
    return "NA"


def strand(raw: dict[str, ty.Any]) -> str:
    for column in STRAND_COLUMNS:
        value = raw.get(column)
        if value and str(value).strip():
            return str(value).strip()
    return "."


def metadata_frame(redi_metadata, genome_build: str) -> pd.DataFrame:
    metadata = pd.read_csv(
        redi_metadata, delimiter="\t", usecols=["Region", "Position", "Ref", "Ed"]
    )

    metadata["region_id"] = metadata.apply(
        lambda row: region_id(row["Region"], row["Position"], genome_build),
        axis=1,
    )
    return metadata


def bed_rows(redi_metadata, genome_build: str) -> ty.Iterable[list[ty.Any]]:
    reader = csv.DictReader(redi_metadata, delimiter="\t")
    for row in reader:
        chrom = row["Region"]
        position = int(row["Position"])
        yield [
            chrom,
            position,
            position + 1,
            region_id(chrom, position, genome_build),
            0,
            strand(row),
            repeat_type(row),
        ]


@attr.s(frozen=True)
class GenomicLocation(object):
    chromosome: str = attr.ib(validator=is_a(str))
    start: int = attr.ib(validator=is_a(int), converter=int)
    stop: int = attr.ib(validator=is_a(int), converter=int)

    @classmethod
    def build(cls, raw):
        chromosome = re.sub("^chr", "", raw.chrom)
        return cls(
            chromosome=chromosome,
            start=raw.start_rel_genome,
            stop=raw.end_rel_genome,
        )


@attr.s(frozen=True)
class RNAEditFeature:
    upi: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int), converter=int)
    repeat_type: str = attr.ib(validator=is_a(str))
    ref: str = attr.ib(validator=is_a(str))
    ed: str = attr.ib(validator=is_a(str))
    start: int = attr.ib(validator=is_a(int), converter=int)
    stop: int = attr.ib(validator=is_a(int), converter=int)
    genomic_location: GenomicLocation = attr.ib(validator=is_a(GenomicLocation))
    genome_build: str = attr.ib(validator=is_a(str))

    @classmethod
    def build(cls, raw_feature):
        upi, taxid = raw_feature.urs_taxid.split("_")
        return cls(
            upi=upi,
            taxid=taxid,
            repeat_type=raw_feature.repeat_type,
            ref=raw_feature.Ref,
            ed=raw_feature.Ed,
            start=raw_feature.start_rel_URS,
            stop=raw_feature.end_rel_URS,
            genomic_location=GenomicLocation.build(raw_feature),
            genome_build=raw_feature.genome_build,
        )

    def writeable(self):
        metadata = attr.asdict(self)

        # metadata['genomic_location']['strand'] = metadata['genomic_location']['strand'].display_int()
        metadata = {
            "repeat_type": self.repeat_type,
            "genomic_location": metadata["genomic_location"],
            "reference": self.ref,
            "edit": self.ed,
            "url": editing_url(
                self.genomic_location.chromosome,
                self.genomic_location.start,
                self.genomic_location.stop,
                self.genome_build,
            ),
        }

        return [
            self.upi,
            self.taxid,
            None,
            self.start,
            self.stop,
            "rna_editing_event",
            json.dumps(metadata),
            "REDIPORTAL",
        ]


def build_bed(redi_metadata, output, genome_build: str) -> None:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    for row in bed_rows(redi_metadata, genome_build):
        writer.writerow(row)


def parse(redi_bedfile, redi_metadata, rnc_bedfile, output, genome_build: str):
    import pybedtools as pbt

    redi_bed = pbt.BedTool(redi_bedfile)
    rnc_bed = pbt.BedTool(rnc_bedfile)

    metadata = metadata_frame(redi_metadata, genome_build)

    ## Make a dataframe from a bed file. Prefix things I want to drop with an underscore, but they all
    ## Need to be present for the conversion to work
    intersection = (
        redi_bed.intersect(rnc_bed, wb=True, s=True)  # , sorted=True
        .to_dataframe(
            names=[
                "chrom",
                "start_rel_genome",
                "end_rel_genome",
                "region_id",
                "score",
                "strand",
                "repeat_type",
                "_rnc_chrom",
                "rnc_exon_start",
                "rnc_exon_end",
                "urs_taxid",
                "urs_length",
                "_rnc_strand",
                "rnc_transcript_start",
                "rnc_transcript_end",
            ],
            index_col=False,
        )
        .drop(
            [
                "_rnc_chrom",
                "_rnc_strand",
            ],
            axis="columns",
        )
    )

    intersection["start_rel_URS"] = np.where(
        intersection["strand"] == "-",
        intersection["rnc_exon_start"]
        - intersection["start_rel_genome"]
        + intersection["urs_length"],
        intersection["start_rel_genome"] - intersection["rnc_transcript_start"],
    )
    intersection["end_rel_URS"] = intersection["start_rel_URS"]
    intersection["genome_build"] = genome_build

    complete_data = intersection.merge(metadata, how="inner", on="region_id")

    def _rows():
        for hit in complete_data.itertuples(index=False):
            yield RNAEditFeature.build(hit).writeable()

    if isinstance(output, (str, Path)):
        path = Path(output)
        if is_parquet():
            with parquet_writer(path, schemas.REDIPORTAL_FEATURES) as writer:
                for row in _rows():
                    writer.writerow(row)
            return
        with path.open("w") as handle:
            writer = csv.writer(handle, delimiter=",")
            for row in _rows():
                writer.writerow(row)
        return

    writer = csv.writer(output, delimiter=",")
    for row in _rows():
        writer.writerow(row)
