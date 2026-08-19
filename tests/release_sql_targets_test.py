# -*- coding: utf-8 -*-

"""
load-data.nf selects a pre/post-release script by globbing
``*__<logical-name>.sql`` for each name that was loaded, where the loaded names
are the keys of ``ENTRY_WRITER_LOAD_TABLES``. Nothing checks that the script's
destination table actually exists -- so mapping a new name pulls in a script
that can die mid-release, after hours of index builds, on "relation ... does not
exist".

That is exactly how 001__karyotypes.sql shipped: it targeted ``ensembl_karyotype``,
a table that has never existed in the RNAcentral schema. It sat unreachable for
years (its parser emits no rows, so the 0-byte csv was dropped by load-data.nf's
``!f.isEmpty()`` filter) until the parquet writer started emitting a non-empty
zero-row file and the name was mapped.

Destination tables are owned by rnacentral-webcode, not this repo, so the
inventory below is a checked-in snapshot verified against the RNAcentral schema.
Adding a destination means confirming the table exists and listing it here.
"""

import re
from pathlib import Path

import pytest

from rnacentral_pipeline import schemas

ROOT = Path(__file__).resolve().parents[1]
RELEASE_DIRS = [
    ROOT / "files" / "import-data" / "pre-release",
    ROOT / "files" / "import-data" / "post-release",
]

# Run unconditionally by load-data.nf rather than selected by a loaded name.
ALWAYS_RUN = {"000__populate_precompute.sql", "999__cleanup.sql"}

# Tables in the RNAcentral schema that release scripts write to. Verified to
# exist; ``rnc_coordinates`` and ``ref_pubmed`` are deliberately absent because
# they do not, which is why the scripts naming them must stay unreachable.
DESTINATION_TABLES = {
    "ensembl_assembly",
    "ensembl_compara",
    "ensembl_coordinate_systems",
    "ensembl_import_tracking",
    "ensembl_pseudogene_exons",
    "ensembl_pseudogene_regions",
    "go_term_annotations",
    "go_term_publication_map",
    "ontology_terms",
    "protein_info",
    "rfam_clans",
    "rfam_go_terms",
    "rfam_models",
    "rnc_accession_sequence_region",
    "rnc_gene_members",
    "rnc_genes",
    "rnc_interactions",
    "rnc_related_sequences",
    "rnc_rna_precomputed",
    "rnc_secondary_structure",
    "rnc_sequence_exons",
    "rnc_sequence_features",
    "rnc_sequence_regions",
    "rnc_taxonomy",
}

WRITE = re.compile(
    r"(?:INSERT\s+INTO|DELETE\s+FROM|TRUNCATE\s+(?:TABLE\s+)?)\s*([\w.]+)"
    r"|(?<!DO )\bUPDATE\s+([\w.]+)\s+SET",
    re.I,
)
CREATE = re.compile(
    r"CREATE\s+(?:TEMP(?:ORARY)?\s+|UNLOGGED\s+)?TABLE\s+(?:IF\s+NOT\s+EXISTS\s+)?([\w.]+)",
    re.I,
)

LOAD_TABLES = {
    m.group(1).lower()
    for m in CREATE.finditer(
        (ROOT / "files" / "schema" / "create_load.sql").read_text()
    )
}


def written_tables(sql: str) -> set:
    """Tables a script writes to, minus the ones it creates for itself."""
    local = {m.group(1).lower().split(".")[-1] for m in CREATE.finditer(sql)}
    found = set()
    for match in WRITE.finditer(sql):
        found.add((match.group(1) or match.group(2)).lower().split(".")[-1])
    return found - local


def reachable_scripts() -> list:
    """The scripts load-data.nf can actually select, by its own glob."""
    stems = {name.replace("_", "-") for name in schemas.ENTRY_WRITER_LOAD_TABLES}
    found = []
    for directory in RELEASE_DIRS:
        for script in sorted(directory.glob("*.sql")):
            _, _, stem = script.stem.partition("__")
            if stem in stems or script.name in ALWAYS_RUN:
                found.append(script)
    return found


REACHABLE = reachable_scripts()


@pytest.mark.parametrize("script", REACHABLE, ids=[s.name for s in REACHABLE])
def test_reachable_release_scripts_target_existing_tables(script):
    known = LOAD_TABLES | DESTINATION_TABLES
    for table in sorted(written_tables(script.read_text())):
        assert table in known, (
            f"{script.name} writes to {table}, which is neither a staging table "
            "in create_load.sql nor a known RNAcentral table"
        )


def test_scripts_naming_missing_tables_stay_unreachable():
    """
    ``rnc_coordinates`` and ``ref_pubmed`` no longer exist. Their scripts are
    harmless only because no loaded name globs them; mapping ``locations`` or
    ``go_publications`` would kill a release the way karyotypes did.
    """
    selectable = {s.name for s in REACHABLE}
    assert "001__locations.sql" not in selectable
    assert "002__go-publications.sql" not in selectable
