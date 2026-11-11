# -*- coding: utf-8 -*-

"""
Checks the text written to the genome-coordinates files from hand-built records, so
they run without a database. The IGV GFF3 is now derived from the GFF3, so these
columns are what both files depend on.
"""

import csv
import io
import json
import logging

import pytest

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import bed, gff3
from rnacentral_pipeline.rnacentral.ftp_export.coordinates import data as coord

EXON_1 = {"exon_start": 100, "exon_stop": 200}
EXON_2 = {"exon_start": 300, "exon_stop": 400}


def raw(**overrides):
    record = {
        "rna_id": "URS1_9606",
        "gene_name": None,
        "description": "a snoRNA",
        "rna_type": "snoRNA",
        "databases": ["ENA", "snoDB"],
        "providing_databases": ["snoDB"],
        "was_mapped": False,
        "identity": None,
        "exons": [EXON_1],
        "assembly_id": "GRCh38",
        "chromosome": "3",
        "strand": 1,
    }
    record.update(overrides)
    return record


def as_csv(*records):
    # The query exports each record as one CSV-quoted JSON field.
    out = io.StringIO()
    writer = csv.writer(out)
    for record in records:
        writer.writerow([json.dumps(record)])
    return io.StringIO(out.getvalue())


def gff(*records):
    out = io.StringIO()
    gff3.from_file(as_csv(*records), out)
    return [line.split("\t") for line in out.getvalue().splitlines()]


def attributes(column):
    return dict(pair.split("=", 1) for pair in column.split(";"))


def test_a_transcript_is_a_transcript_line_then_one_line_per_exon():
    rows = gff(raw(exons=[EXON_1, EXON_2]))

    assert rows[0] == ["##gff-version 3"]
    transcript, exon_1, exon_2 = rows[1:]
    assert transcript[:8] == [
        "3",
        "RNAcentral",
        "transcript",
        "100",
        "400",
        ".",
        "+",
        ".",
    ]
    assert [exon_1[2:5], exon_2[2:5]] == [
        ["noncoding_exon", "100", "200"],
        ["noncoding_exon", "300", "400"],
    ]
    assert attributes(transcript[8])["ID"] == "URS1_9606.0"
    assert attributes(exon_2[8])["ID"] == "URS1_9606.0:ncRNA_exon2"
    assert attributes(exon_2[8])["Parent"] == "URS1_9606.0"


def test_an_expert_database_region_names_its_providers():
    transcript = gff(raw())[1]

    attrs = attributes(transcript[8])
    assert attrs["source"] == "expert-database"
    assert attrs["providing_databases"] == "snoDB"
    assert attrs["databases"] == "ENA,snoDB"


def test_a_mapped_region_reports_its_identity_to_two_decimal_places():
    transcript = gff(raw(was_mapped=True, identity=0.987, providing_databases=[]))[1]

    attrs = attributes(transcript[8])
    assert attrs["source"] == "alignment"
    assert attrs["identity"] == "0.99"
    assert "providing_databases" not in attrs


def test_a_transcript_points_at_its_gene():
    transcript = gff(raw(gene_name="G1"))[1]

    assert attributes(transcript[8])["Parent"] == "G1"


def test_a_gene_uses_its_own_start_and_end():
    gene = raw(
        rna_id=None,
        gene_name="GENE1",
        rna_type="lncRNA",
        providing_databases=[],
        exons=[],
        gene_start=10,
        gene_end=50,
    )

    row = gff(gene)[1]

    assert row[2:7] == ["predicted_gene", "10", "50", ".", "+"]
    assert attributes(row[8])["ID"] == "GENE1"


def test_a_reverse_strand_is_written_as_a_minus():
    assert gff(raw(strand=-1))[1][6] == "-"


def test_regions_without_a_length_are_skipped(caplog):
    flat = raw(exons=[{"exon_start": 100, "exon_stop": 100}])

    with caplog.at_level(logging.WARNING):
        rows = gff(flat, raw(rna_id="URS2_9606"))

    assert [r[2] for r in rows[1:]] == ["transcript", "noncoding_exon"]
    assert attributes(rows[1][8])["Name"] == "URS2_9606"
    assert "Skipping region" in caplog.text


def test_writing_no_features_is_an_error_unless_allowed():
    with pytest.raises(ValueError):
        gff3.write_gff_text(iter([]), io.StringIO())

    assert gff3.write_gff_text(iter([]), io.StringIO(), allow_no_features=True) is False


def test_the_gff3_header_can_be_left_off_for_a_second_file():
    features = gff3.regions_as_features(coord.from_file(as_csv(raw())))
    out = io.StringIO()

    assert gff3.write_gff_text(features, out, header=False) is True
    assert not out.getvalue().startswith("##gff-version")


def bed_rows(function, *records):
    out = io.StringIO()
    function(as_csv(*records), out)
    return [line.split("\t") for line in out.getvalue().splitlines()]


def test_bed_is_zero_based_with_a_block_per_exon():
    (row,) = bed_rows(bed.from_json, raw(exons=[EXON_1, EXON_2]))

    assert row == [
        "chr3", "99", "400", "URS1_9606", "0", "+", "99", "400",
        "63,125,151", "2", "101,101", "0,200", ".", "snoRNA", "ENA,snoDB",
    ]  # fmt: skip


@pytest.mark.parametrize(
    "name,expected", [("MT", "chrM"), ("chrMT", "chrM"), ("X", "chrX")]
)
def test_bed_chromosome_names(name, expected):
    (row,) = bed_rows(bed.from_json, raw(chromosome=name))

    assert row[0] == expected


def test_exploded_bed_has_one_row_per_exon_with_the_transcript_length():
    rows = bed_rows(bed.from_json_exploded, raw(exons=[EXON_1, EXON_2]))

    assert rows == [
        ["chr3", "99", "200", "URS1_9606", "202", "+", "99", "400"],
        ["chr3", "299", "400", "URS1_9606", "202", "+", "99", "400"],
    ]


def test_bed_refuses_overlapping_exons():
    overlapping = raw(exons=[EXON_1, {"exon_start": 100, "exon_stop": 150}])

    with pytest.raises(AssertionError):
        bed_rows(bed.from_json, overlapping)


def test_genes_are_left_out_of_the_bed_file():
    gene = raw(
        rna_id=None,
        gene_name="G",
        providing_databases=[],
        exons=[],
        gene_start=1,
        gene_end=5,
    )

    assert len(bed_rows(bed.from_json, gene, raw())) == 1
