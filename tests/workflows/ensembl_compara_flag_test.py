from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def text(path):
    return (ROOT / path).read_text()


# Compara files are not always published; a compara-only run loads them once they are.
def test_compara_is_off_by_default_with_its_own_flag():
    config = text("config/databases.config")
    compara = config.split("compara {", 1)[1].split("}", 1)[0]
    assert "run = false" in compara


def test_compara_flag_alone_fetches_compara_and_parses_ensembl_metadata():
    fetch = text("workflows/metadata/ensembl.nf").split("process fetch_compara {", 1)[1]
    assert "params.databases.ensembl?.compara?.run" in fetch.split("script:", 1)[0]
    assert "params.databases.ensembl?.compara?.run" in text(
        "workflows/parse-metadata.nf"
    )


def test_compara_flag_alone_triggers_the_release_that_loads_it():
    utils = text("lib/Utils.groovy")
    should_release = utils.split("static def should_release", 1)[1].split("}\n", 1)[0]
    assert "will_run(compara)" in should_release


def test_compara_flag_does_not_start_the_species_parse_or_proteins():
    ensembl = text("workflows/metadata/ensembl.nf")
    assemblies = ensembl.split("process assemblies {", 1)[1].split("script:", 1)[0]
    proteins = ensembl.split("process proteins {", 1)[1].split("script:", 1)[0]
    assert "compara" not in assemblies and "compara" not in proteins
