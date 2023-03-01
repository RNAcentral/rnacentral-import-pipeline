import csv
import json
import re

import attr
import pandas as pd
import pybedtools as pbt
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases import data

REDI_BASE_URL = "http://srv00.recas.ba.infn.it/cgi/atlas/getpage_dev.py?query1=chr{0}:{1}-{2}&query10=hg38&query9=hg"


@attr.s(frozen=True)
class GenomicLocation(object):
    chromosome = attr.ib(validator=is_a(str))
    start = attr.ib(validator=is_a(int), converter=int)
    stop = attr.ib(validator=is_a(int), converter=int)

    @classmethod
    def build(cls, raw):
        chromosome = re.sub("^chr", "", raw.chrom)
        return cls(
            chromosome=chromosome,
            start=raw.start_rel_genome,
            stop=raw.end_rel_genome,
        )


@attr.s(frozen=True)
class RNAEditFeature(object):
    upi = attr.ib(validator=is_a(str))
    taxid = attr.ib(validator=is_a(int), converter=int)
    repeat_type = attr.ib(validator=is_a(str))
    ref = attr.ib(validator=is_a(str))
    ed = attr.ib(validator=is_a(str))
    start = attr.ib(validator=is_a(int), converter=int)
    stop = attr.ib(validator=is_a(int), converter=int)
    genomic_location = attr.ib(validator=is_a(GenomicLocation))

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
        )

    def writeable(self):
        metadata = attr.asdict(self)
        metadata = {
            "repeat_type": self.repeat_type,
            "genomic_location": metadata["genomic_location"],
            "reference": self.ref,
            "edit": self.ed,
            "url": REDI_BASE_URL.format(
                self.genomic_location.chromosome,
                self.genomic_location.start,
                self.genomic_location.stop,
            ),
        }

        return [
            self.upi,
            self.taxid,
            None,
            self.start,
            self.stop,
            "RNA_editing_event",
            json.dumps(metadata),
            "REDIPORTAL",
        ]


def parse(redi_bedfile, redi_metadata, rnc_bedfile, output):
    redi_bed = pbt.BedTool(redi_bedfile)
    rnc_bed = pbt.BedTool(rnc_bedfile)

    metadata = pd.read_csv(
        redi_metadata, delimiter="\t", usecols=["Region", "Position", "Ref", "Ed"]
    )

    metadata["region_id"] = (
        metadata["Region"] + "_" + metadata["Position"].astype(str) + "_hg38"
    )

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
                "_rnc_score",
                "_rnc_strand",
                "rnc_transcript_start",
                "rnc_transcript_end",
            ],
            index_col=False,
        )
        .drop(
            [
                "_rnc_chrom",
                "_rnc_score",
                "_rnc_strand",
            ],
            axis="columns",
        )
    )

    intersection["start_rel_URS"] = (
        intersection["start_rel_genome"] - intersection["rnc_transcript_start"]
    )
    intersection["end_rel_URS"] = (
        intersection["end_rel_genome"] - intersection["rnc_transcript_start"]
    )

    complete_data = intersection.merge(metadata, how="inner", on="region_id")

    writer = csv.writer(output, delimiter=",")
    for hit in complete_data.itertuples(index=False):
        ef = RNAEditFeature.build(hit)
        writer.writerow(ef.writeable())
