import io
import json
from pathlib import Path

from rnacentral_pipeline.databases.rediportal import parser

FIXTURE = Path(__file__).with_name("TABLE1_hg38_v3_subset.tsv")


def test_build_bed_uses_region_position_and_genome_build():
    metadata = FIXTURE.open()
    output = io.StringIO()

    try:
        parser.build_bed(metadata, output, "hg38")
    finally:
        metadata.close()

    assert output.getvalue() == "chr1\t87158\t87159\tchr1_87158_hg38\t0\t-\tSINE/AluJo\n"


def test_build_bed_falls_back_when_repeat_columns_are_missing():
    metadata = io.StringIO(
        "Region\tPosition\tRef\tEd\tStrand\n"
        "chr2\t99\tA\tI\t-\n"
    )
    output = io.StringIO()

    parser.build_bed(metadata, output, "mm10")

    assert output.getvalue() == "chr2\t99\t100\tchr2_99_mm10\t0\t-\tNA\n"


def test_feature_metadata_uses_requested_genome_build():
    feature = parser.RNAEditFeature(
        upi="URS0001",
        taxid=9606,
        repeat_type="ALU",
        ref="A",
        ed="G",
        start=5,
        stop=5,
        genomic_location=parser.GenomicLocation("1", 42, 43),
        genome_build="mm10",
    )

    writeable = feature.writeable()
    metadata = json.loads(writeable[6])

    assert (
        metadata["url"]
        == "https://rediportal.cloud.ba.infn.it/cgi/atlas/getpage_dev.py?"
        "query1=chr1%3A42-43&query10=mm10&query9=mm"
    )
