# -*- coding: utf-8 -*-

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "split_by_species", ROOT / "bin" / "split_by_species.py"
)
split_by_species = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(split_by_species)

SPECIES_CSV = (
    "9606,Homo sapiens,cellular organisms;Eukaryota;Metazoa;Vertebrata\n"
    "10090,Mus musculus,cellular organisms;Eukaryota;Metazoa;Vertebrata\n"
    "562,Escherichia coli,cellular organisms;Bacteria\n"
)


@pytest.fixture
def species_csv(tmp_path):
    path = tmp_path / "species.csv"
    path.write_text(SPECIES_CSV)
    return path


def read(path):
    return Path(path).read_text()


@pytest.mark.utils
def test_records_are_written_to_a_file_per_species_under_their_kingdom(
    tmp_path, species_csv
):
    fasta = tmp_path / "all.fasta"
    fasta.write_text(
        ">URS0001_9606 human\nACGU\n>URS0002_9606 human\nGGGG\n>URS0003_562 bug\nUUUU\n"
    )

    split_by_species.split(species_csv, fasta, tmp_path / "out")

    assert (
        read(tmp_path / "out/vertebrates/homo_sapiens.fasta")
        == ">URS0001_9606 human\nACGU\n>URS0002_9606 human\nGGGG\n"
    )
    assert (
        read(tmp_path / "out/bacteria/escherichia_coli.fasta")
        == ">URS0003_562 bug\nUUUU\n"
    )


@pytest.mark.utils
def test_taxa_that_are_not_listed_are_skipped(tmp_path, species_csv):
    fasta = tmp_path / "all.fasta"
    fasta.write_text(">URS0001_1 unlisted\nACGU\n>URS0002_9606 human\nGGGG\n")

    split_by_species.split(species_csv, fasta, tmp_path / "out")

    assert (
        read(tmp_path / "out/vertebrates/homo_sapiens.fasta")
        == ">URS0002_9606 human\nGGGG\n"
    )


@pytest.mark.utils
def test_a_taxid_that_comes_back_later_is_an_error_not_a_truncated_file(
    tmp_path, species_csv
):
    # The input must be grouped by taxid; reopening a species file would overwrite it.
    fasta = tmp_path / "all.fasta"
    fasta.write_text(
        ">URS0001_9606 a\nACGU\n>URS0003_562 b\nUUUU\n>URS0002_9606 c\nGGGG\n"
    )

    with pytest.raises(SystemExit):
        split_by_species.split(species_csv, fasta, tmp_path / "out")


@pytest.mark.utils
def test_species_sharing_a_name_are_told_apart_by_taxid(tmp_path):
    csv = tmp_path / "species.csv"
    csv.write_text("1,Same name,Bacteria\n2,Same name,Bacteria\n")
    fasta = tmp_path / "all.fasta"
    fasta.write_text(">URS0001_1 a\nA\n>URS0002_2 b\nC\n")

    split_by_species.split(csv, fasta, tmp_path / "out")

    assert (tmp_path / "out/bacteria/same_name_1.fasta").exists()
    assert (tmp_path / "out/bacteria/same_name_2.fasta").exists()
