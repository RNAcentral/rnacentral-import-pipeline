# -*- coding: utf-8 -*-

import io

from rnacentral_pipeline.databases.ensembl.urls import GenesetUrls
from rnacentral_pipeline.rnacentral.genome_mapping import remote_files


def geneset(gca, genome_url):
    return GenesetUrls(
        species="Danio_rerio",
        taxid=7955,
        embl_url="B/genes.embl.gz",
        gff_url=f"B/{gca}/genes.gff3.gz",
        accession=gca,
        genome_url=genome_url,
    )


def test_resolves_genome_and_gff3_by_gca():
    selected = [geneset("GCA_1.1", "B/GCA_1.1/unmasked.fa.bgz")]
    assemblies = io.StringIO("danio_rerio,GRCz11,7955,GCA_1.1\n")
    assert list(remote_files.rows(selected, assemblies)) == [
        [
            "danio_rerio",
            "GRCz11",
            "7955",
            "B/GCA_1.1/unmasked.fa.bgz",
            "B/GCA_1.1/genes.gff3.gz",
        ]
    ]


# Assemblies Ensembl no longer serves, or that list no genome, are skipped, not failed.
def test_skips_assemblies_without_a_genome_in_species_json():
    selected = [geneset("GCA_1.1", None)]
    assemblies = io.StringIO(
        "danio_rerio,GRCz11,7955,GCA_1.1\nold_species,v1,1,GCA_9.1\n"
    )
    assert list(remote_files.rows(selected, assemblies)) == []
