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


class FakeFtp:
    def __init__(self, names):
        self.names = names

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def login(self):
        return None

    def nlst(self, _remote_dir):
        return self.names


def test_resolve_ftp_urls_resolves_wildcard_matches(monkeypatch):
    monkeypatch.setattr(
        url_helpers,
        "FTP",
        lambda _host: FakeFtp(
            [
                "/pub/release-115/embl/ursus_americanus/Ursus_americanus.ASM334442v1.115.nonchromosomal.dat.gz",
            ]
        ),
    )

    found = url_helpers.resolve_ftp_urls(
        "ftp://ftp.ensembl.org/pub/release-115/embl/ursus_americanus/Ursus_americanus.ASM334442v1.115.*.dat.gz"
    )

    assert found == [
        "ftp://ftp.ensembl.org/pub/release-115/embl/ursus_americanus/Ursus_americanus.ASM334442v1.115.nonchromosomal.dat.gz",
    ]


def test_resolve_ftp_urls_falls_back_to_lowercase_assembly(monkeypatch):
    monkeypatch.setattr(
        url_helpers,
        "FTP",
        lambda _host: FakeFtp(
            [
                "/pub/current/plants/embl/arabidopsis_thaliana/Arabidopsis_thaliana.tair10.62.chromosome.1.dat.gz",
            ]
        ),
    )

    found = url_helpers.resolve_ftp_urls(
        "ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.*.dat.gz"
    )

    assert found == [
        "ftp://ftp.ensemblgenomes.org/pub/current/plants/embl/arabidopsis_thaliana/Arabidopsis_thaliana.tair10.62.chromosome.1.dat.gz",
    ]


def test_resolve_ftp_urls_falls_back_to_single_nonchromosomal_file(monkeypatch):
    monkeypatch.setattr(
        url_helpers,
        "FTP",
        lambda _host: FakeFtp(
            [
                "/pub/release-115/embl/sus_scrofa_gca018555405v1/"
                "Sus_scrofa_gca018555405v1.ASM1855540v1.115.nonchromosomal.dat.gz",
            ]
        ),
    )

    found = url_helpers.resolve_ftp_urls(
        "ftp://ftp.ensembl.org/pub/release-115/embl/sus_scrofa_gca018555405v1/"
        "Sus_scrofa_nihs_2020.ASM1855540v1.115.*.dat.gz"
    )

    assert found == [
        "ftp://ftp.ensembl.org/pub/release-115/embl/sus_scrofa_gca018555405v1/"
        "Sus_scrofa_gca018555405v1.ASM1855540v1.115.nonchromosomal.dat.gz",
    ]


def test_resolve_ftp_urls_does_not_guess_when_multiple_nonchromosomal_files_exist(monkeypatch):
    monkeypatch.setattr(
        url_helpers,
        "FTP",
        lambda _host: FakeFtp(
            [
                "/pub/release-115/embl/example/one.nonchromosomal.dat.gz",
                "/pub/release-115/embl/example/two.nonchromosomal.dat.gz",
            ]
        ),
    )

    found = url_helpers.resolve_ftp_urls(
        "ftp://ftp.ensembl.org/pub/release-115/embl/example/missing.*.dat.gz"
    )

    assert found == []
