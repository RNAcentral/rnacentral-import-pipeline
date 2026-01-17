# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

"""
Parser for CIRCpedia V3 circular RNA data.

CIRCpedia V3 provides comprehensive circular RNA annotations from
multiple species with expression data and genomic coordinates.
"""

import logging
import typing as ty
from pathlib import Path

import polars as pl
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.circpedia import helpers
from rnacentral_pipeline.databases.helpers import phylogeny as phy

LOGGER = logging.getLogger(__name__)

# RNA type for circular RNA (SO:0000593)
CIRCULAR_RNA_SO_TERM = "SO:0000593"


def parse_csv_row(
    row: pl.Series,
    taxonomy: SqliteDict,
    assembly_id: ty.Optional[str] = None,
) -> ty.Optional[Entry]:
    """
    Parse a single row from CIRCpedia CSV data.

    Args:
        row: Polars Series representing one circRNA
        taxonomy: SqliteDict with taxonomy data
        assembly_id: Genome assembly ID (e.g., GRCh38)

    Returns:
        Entry object or None if parsing fails
    """
    try:
        # Get circRNA ID (required)
        if "circid" not in row:
            LOGGER.warning("Missing circid field in row")
            return None

        circ_id = row["circid"]
        if circ_id is None or (isinstance(circ_id, float) and circ_id != circ_id):
            LOGGER.warning("Invalid circid in row")
            return None

        circ_id = str(circ_id).strip()

        # Get species and taxonomy ID (required)
        if "species" not in row:
            LOGGER.warning(f"Missing species field for {circ_id}")
            return None

        species_str = str(row["species"]).strip() if row["species"] else None
        if not species_str:
            LOGGER.warning(f"Empty species for {circ_id}")
            return None

        taxid = helpers.get_ncbi_taxid(taxonomy, species_str)
        if not taxid:
            LOGGER.warning(f"Could not determine taxid for species: {species_str}")
            return None

        # Parse genomic location (required for circRNAs)
        if "location" not in row or not row["location"]:
            LOGGER.warning(f"Missing location for {circ_id}")
            return None

        location = helpers.parse_genomic_location(row["location"])
        if not location:
            LOGGER.warning(f"Could not parse location for {circ_id}: {row['location']}")
            return None

        # Get strand information
        strand = 0
        if "strand" in row and row["strand"]:
            strand = helpers.parse_strand(row["strand"])

        # Parse exon positions if available
        exons_list = []
        if "exon_positions" in row and row["exon_positions"]:
            exons_list = helpers.parse_exon_positions(row["exon_positions"])
        elif "exonstart_exonend" in row and row["exonstart_exonend"]:
            # Alternative column name
            exons_list = helpers.parse_exon_positions(row["exonstart_exonend"])

        # Get sequence if provided, otherwise use placeholder
        # Circular RNA sequences may need to be constructed from exons
        sequence = ""
        if "sequence" in row and row["sequence"]:
            sequence = str(row["sequence"]).upper().replace('U', 'T')

        # If no sequence provided, create a placeholder
        # In production, you would fetch from genome or database
        if not sequence:
            # Use 'N' * length as placeholder - this should be replaced
            # with actual sequence extraction in production
            seq_length = location["end"] - location["start"] + 1
            sequence = "N" * min(seq_length, 100)  # Limit placeholder size
            LOGGER.warning(
                f"No sequence provided for {circ_id}, using placeholder. "
                "Sequence should be extracted from genome."
            )

        # Get gene name
        gene_name = helpers.gene(row)

        # Build the Entry
        entry = Entry(
            primary_id=helpers.primary_id(circ_id),
            accession=helpers.accession(circ_id, location),
            ncbi_tax_id=taxid,
            database="CIRCPEDIA",
            sequence=sequence,
            regions=helpers.regions(location, strand, exons_list, assembly_id),
            rna_type=CIRCULAR_RNA_SO_TERM,
            url=helpers.url(circ_id),
            seq_version=helpers.seq_version(row),
            note_data=helpers.note_data(row, location),
            chromosome=helpers.chromosome(location),
            species=phy.species(taxid, taxonomy),
            common_name=phy.common_name(taxid, taxonomy),
            lineage=phy.lineage(taxid, taxonomy),
            gene=gene_name,
            product=helpers.product(row, gene_name),
            description=helpers.description(taxonomy, taxid, gene_name),
            mol_type="genomic DNA",
            references=helpers.references(),
        )

        return entry

    except phy.FailedTaxonId as e:
        LOGGER.warning(f"Failed to get taxonomy info: {e}")
        return None
    except phy.UnknownTaxonId as e:
        LOGGER.warning(f"Unknown taxonomy ID: {e}")
        return None
    except Exception as e:
        LOGGER.error(f"Error parsing row: {e}", exc_info=True)
        return None


def parse(
    csv_file: ty.Union[str, Path],
    taxonomy_file: Path,
    assembly_id: ty.Optional[str] = None,
) -> ty.Iterable[Entry]:
    """
    Parse CIRCpedia V3 CSV file.

    Expected CSV columns (case-insensitive):
        - circid (required): CIRCpedia circular RNA ID
        - species (required): Species name (e.g., "Homo sapiens")
        - location (required): Genomic location (e.g., "chr1:12345-67890")
        - strand (optional): Strand ('+' or '-')
        - gene (optional): Host gene name
        - sequence (optional): RNA sequence
        - exon_positions or exonstart_exonend (optional): Exon coordinates
        - fpm (optional): Expression level (Fragments Per Million)
        - cell_line (optional): Cell line or tissue
        - sequencing_type (optional): Type of sequencing
        - conservation (optional): Conservation information
        - rnase_r_enrichment (optional): RNase R enrichment fold change

    Args:
        csv_file: Path to CIRCpedia CSV file
        taxonomy_file: Path to taxonomy database file
        assembly_id: Genome assembly ID (e.g., GRCh38, GRCm39)

    Yields:
        Entry objects for each circular RNA
    """
    LOGGER.info(f"Parsing CIRCpedia data from {csv_file}")

    # Load taxonomy data
    taxonomy = SqliteDict(str(taxonomy_file))

    # Read CSV with polars for efficient processing
    try:
        # Read CSV and normalize column names to lowercase
        df = pl.read_csv(
            csv_file,
            infer_schema_length=10000,
            ignore_errors=True,
        )

        # Normalize column names to lowercase for easier matching
        df = df.rename({col: col.lower() for col in df.columns})

        LOGGER.info(f"Read {len(df)} rows from CSV")
        LOGGER.info(f"Columns: {df.columns}")

        # Check for required columns
        required_cols = ["circid", "species", "location"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Process each row
        parsed_count = 0
        error_count = 0

        for row in df.iter_rows(named=True):
            entry = parse_csv_row(row, taxonomy, assembly_id)
            if entry:
                yield entry
                parsed_count += 1
            else:
                error_count += 1

            # Log progress periodically
            if (parsed_count + error_count) % 10000 == 0:
                LOGGER.info(
                    f"Processed {parsed_count + error_count} rows: "
                    f"{parsed_count} parsed, {error_count} errors"
                )

        LOGGER.info(
            f"Parsing complete: {parsed_count} entries created, "
            f"{error_count} rows failed"
        )

    except Exception as e:
        LOGGER.error(f"Error reading CSV file: {e}", exc_info=True)
        raise
    finally:
        taxonomy.close()
