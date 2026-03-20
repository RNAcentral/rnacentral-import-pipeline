# -*- coding: utf-8 -*-

from rnacentral_pipeline.databases.ensembl import url_helpers


def test_lowercase_assembly_in_url_changes_only_the_assembly_component():
    url = (
        "ftp://ftp.ensemblgenomes.org/pub/current/plants/gff3/"
        "arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz"
    )

    found = url_helpers.lowercase_assembly_in_url(url)

    assert found == (
        "ftp://ftp.ensemblgenomes.org/pub/current/plants/gff3/"
        "arabidopsis_thaliana/Arabidopsis_thaliana.tair10.62.gff3.gz"
    )


def test_lowercase_assembly_in_url_preserves_suffixes_after_release():
    url = (
        "ftp://ftp.ensembl.org/pub/release-110/embl/"
        "homo_sapiens/Homo_sapiens.GRCh38.110.1.dat.gz"
    )

    found = url_helpers.lowercase_assembly_in_url(url)

    assert found == (
        "ftp://ftp.ensembl.org/pub/release-110/embl/"
        "homo_sapiens/Homo_sapiens.grch38.110.1.dat.gz"
    )
