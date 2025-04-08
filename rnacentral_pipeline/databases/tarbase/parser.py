import logging
import os
import sys
import typing as ty

import polars as pl
import psycopg2 as pg
from psycopg2.extras import execute_batch

from rnacentral_pipeline.databases.data.entry import Entry
from rnacentral_pipeline.databases.data.related import RelatedEvidence, RelatedSequence
from rnacentral_pipeline.databases.helpers import phylogeny as phy

logger = logging.getLogger(__name__)

from rnacentral_pipeline.rnacentral.notify.slack import send_notification


def build_entry(row: ty.Dict[str, str]) -> Entry:
    primary_id = row["mirna_name"]
    ## Use this accession form to match existing accessions
    accession = f"TARBASE:{row['mirna_name']}"

    ## This takes you to the specific search results for the interation on this row
    url = f"https://dianalab.e-ce.uth.gr/tarbasev9/interactions?gene={row['gene_name'][0]}&mirna={row['mirna_name']}"

    related_sequences = []
    for (
        related_gene,
        exp_meth,
    ) in zip(row["gene_id"], row["experimental_method"]):

        evidence = RelatedEvidence(methods=[exp_meth])

        related_sequences.append(
            RelatedSequence(
                sequence_id=f"ENSEMBL:{related_gene}",
                relationship="target_protein",
                evidence=evidence,
            )
        )

    ## Make this as minimal as possible, don't fill everything
    ent = Entry(
        primary_id=primary_id,
        accession=accession,
        ncbi_tax_id=int(row["taxid"]),
        database="TARBASE",
        sequence=row["sequence"].replace("U", "T"),
        regions=[],
        rna_type="SO:0000276",  # miRNA - tarBase is all miRNA
        url=url,
        seq_version=1,
        optional_id=row["mirna_id"],
        description=f"{row['species']} ({phy.common_name(int(row['taxid'])) }) {row['mirna_name']}",
        note_data={"url": url},
        xref_data={},
        related_sequences=related_sequences,
    )

    return ent


def parse(filepath: str) -> ty.List[Entry]:
    """
    Parse the provided tsv file into entry objects we can put in the database

    Some notes on the TSV:
    - mirna_name and mirna_id are the miRNA identifiers we will hold and should look for
    - gene_name, gene_id, transcript_name and transcript_id all reter to the protein coding gene
        targeted by the miRNA

    """
    data = pl.read_csv(filepath, separator="\t", null_values=["NA"]).unique()

    ## The taxa query will be null for the viruses provided by tarbase, so we make some substitutions
    ## These are the NCBI names based on searching NBCI for the given species, or looking in RNAc
    ## for interacting miRNAs
    data = data.with_columns(
        pl.col("species").replace(
            ["EBV", "HIV-1", "KSHV", "MHV68"],
            [
                "human gammaherpesvirus 4",
                "Human immunodeficiency virus 1",
                "Human gammaherpesvirus 8",
                "Murid gammaherpesvirus 4",
            ],
        )
    )

    species_mentioned = tuple(data.get_column("species").unique().to_list())
    ## Try to catch when we miss something - raise a warning
    n_species = len(species_mentioned)

    taxa_query = """
    SELECT name as species, id as taxid FROM rnc_taxonomy
    WHERE name IN %(names)s
    AND replaced_by IS NULL
    """

    ## Have to use a psycopg connection to allow it to mogrify names into the query
    conn = pg.connect(os.getenv("PGDATABASE"))
    cur = conn.cursor()
    taxa_results = pl.read_database(
        query=taxa_query,
        connection=conn,
        execute_options={"parameters": {"names": species_mentioned}},
    )
    if taxa_results.height < n_species:
        if taxa_results.height > 0:
            missing = (
                data.join(taxa_results, on="species", how="anti")
                .get_column("species")
                .unique()
                .to_list()
            )
        else:
            missing = species_mentioned
        message = (
            "Tarbase taxonomy query returned fewer results than we asked for species\n"
            "We are probably missing some entries! Ignoring for now.\n"
            f"Missing: {', '.join(missing)}"
        )
        logger.warning(message)

        ## Only warn on slack if not testing
        if not "pytest" in sys.modules:
            send_notification("Tarbase parser warning", message)

    data = data.join(taxa_results, on="species", how="left").drop_nulls("taxid")

    grouped_data = (
        data.group_by("mirna_id", maintain_order=True)
        .agg(
            pl.col("mirna_name").first(),
            pl.col("gene_id"),
            pl.col("gene_name"),
            pl.col("experimental_method"),
            pl.col("article_pubmed_id"),
            pl.col("taxid").first(),
            pl.col("species").first(),
        )
        .with_columns(
            pl.col("gene_id").list.unique(maintain_order=True),
            pl.col("gene_name").list.unique(maintain_order=True),
            pl.col("experimental_method").list.unique(maintain_order=True),
            pl.col("article_pubmed_id").list.unique(maintain_order=True),
        )
    )

    accessions = [
        (f"MIRBASE:{mirna_id}",)
        for mirna_id in grouped_data.get_column("mirna_id").to_list()
    ]
    cur.execute("CREATE TEMPORARY TABLE temp_accessions_tarbase (ac TEXT)")
    execute_batch(cur, "INSERT INTO temp_accessions_tarbase VALUES (%s)", accessions)

    sequence_query = """
    WITH temp_xref_data AS (
        SELECT upi, xref.ac FROM
        xref JOIN temp_accessions_tarbase ON xref.ac = temp_accessions_tarbase.ac
        WHERE deleted = 'N'
    )
    SELECT
        SPLIT_PART(ac, ':', 2) AS mirna_id,
        COALESCE(seq_short, seq_long) AS sequence
    FROM rna
    JOIN temp_xref_data ON rna.upi = temp_xref_data.upi
    """

    sequence_results = pl.read_database(
        query=sequence_query,
        connection=conn,
        execute_options={
            "parameters": {
                "accessions": accessions,
            }
        },
    )

    grouped_data = grouped_data.join(sequence_results, on="mirna_id").unique()

    entries = [build_entry(r) for r in grouped_data.iter_rows(named=True)]

    return entries
