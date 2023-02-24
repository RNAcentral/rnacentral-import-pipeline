import attr
import pandas as pd
import pybedtools as pbt
from attr.validators import instance_of as is_a


@attr.s(frozen=True)
class GenomicLocation(object):
    chromosome = attr.ib(validator=is_a(str))
    start = attr.ib(validator=is_a(int), converter=int)
    stop = attr.ib(validator=is_a(int), converter=int)

    @classmethod
    def build(cls, raw):
        chromosome = re.sub("^chr", "", raw["chrom"])
        return cls(
            chromosome=chromosome,
            start=raw["start"],
            stop=raw["end"],
        )


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
        redi_bed.intersect(rnc_bed, wb=True, s=True, sorted=True)
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
                "rnc_start",
                "rnc_end",
                "urs_taxid",
                "_rnc_score",
                "_rnc_strand",
                "_rnc_start",
                "_rnc_end",
                "_rnc_bed_rgb",
                "rnc_num_exons",
                "rnc_exon_sizes",
                "rnc_exon_starts",
                "_dot",
                "_rnc_rna_type",
                "_rnc_providing_dbs",
            ],
            index_col=False,
        )
        .drop(
            [
                "_rnc_chrom",
                "_rnc_score",
                "_rnc_strand",
                "_rnc_start",
                "_rnc_end",
                "_rnc_bed_rgb",
                "_dot",
                "_rnc_rna_type",
                "_rnc_providing_dbs",
            ],
            axis="columns",
        )
    )

    # intersection["start_rel_URS"] = (
    #     intersection["start_rel_genome"] - intersection["rnc_start"]
    # )
    # print(intersection.head())
    # print(metadata.head())
    # complete_data = intersection.merge(metadata, how="inner", on="region_id")

    # complete_data["rnc_exon_starts"] = complete_data["rnc_exon_starts"].apply(lambda x: [int(a) for a in x.split(',')])
    # complete_data["rnc_exon_sizes"] = complete_data["rnc_exon_sizes"].apply(lambda x: [int(a) for a in x.split(',')])
    # complete_data = complete_data.explode(["rnc_exon_starts","rnc_exon_sizes"])
    # complete_data["exon_start_rel_genome"] = (complete_data["rnc_exon_starts"] + complete_data["rnc_start"].astype(int)).astype(int)
    # complete_data["exon_end_rel_genome"] = (complete_data["exon_start_rel_genome"] + complete_data["rnc_exon_sizes"]).astype(int)
    # print(complete_data.dtypes)

    # hits = complete_data.loc[complete_data[["start_rel_genome","exon_start_rel_genome"]].max().lt(complete_data[["end_rel_genome","exon_end_rel_genome"]].min(axis=1))]
    # print(hits.head())
