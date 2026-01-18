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

Data format:
- Annotation file: Tab-delimited text (.txt)
- Sequence file: FASTA format (.fa)
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


def load_fasta_sequences(fasta_file: ty.Union[str, Path]) -> ty.Dict[str, str]:
    """
    Load sequences from FASTA file into a dictionary.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        Dictionary mapping circRNA ID to sequence
    """
    sequences = {}
    current_id = None
    current_seq = []

    LOGGER.info(f"Loading sequences from {fasta_file}")

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id and current_seq:
                    sequences[current_id] = ''.join(current_seq)
                # Start new sequence
                current_id = line[1:].strip()  # Remove '>'
                current_seq = []
            elif current_id:
                current_seq.append(line)

        # Save last sequence
        if current_id and current_seq:
            sequences[current_id] = ''.join(current_seq)

    LOGGER.info(f"Loaded {len(sequences)} sequences from FASTA")
    return sequences


def parse_tsv_row(
    row: ty.Dict[str, ty.Any],
    sequences: ty.Dict[str, str],
    taxonomy: SqliteDict,
    assembly_id: ty.Optional[str] = None,
) -> ty.Optional[Entry]:
    """
    Parse a single row from CIRCpedia TSV data.

    Args:
        row: Dictionary representing one circRNA annotation
        sequences: Dictionary mapping circRNA ID to sequence
        taxonomy: SqliteDict with taxonomy data
        assembly_id: Genome assembly ID (e.g., GRCh38)

    Returns:
        Entry object or None if parsing fails
    """
    try:
        # Get circRNA ID (required)
        circ_id = row.get('circID') or row.get('circid')
        if not circ_id or (isinstance(circ_id, float) and circ_id != circ_id):
            LOGGER.warning("Missing or invalid circID in row")
            return None

        circ_id = str(circ_id).strip()

        # Get species and taxonomy ID (required)
        species_str = row.get('species')
        if not species_str:
            LOGGER.warning(f"Missing species field for {circ_id}")
            return None

        species_str = str(species_str).strip()
        # Replace hyphens with spaces for taxonomy lookup
        species_lookup = species_str.replace('-', ' ')

        taxid = helpers.get_ncbi_taxid(taxonomy, species_lookup)
        if not taxid:
            LOGGER.warning(f"Could not determine taxid for species: {species_str}")
            return None

        # Parse genomic location from combined Location field
        # Format: "V:15874634-15876408(-)" or "chr1:100-200(+)"
        location_str = row.get('Location') or row.get('location')
        if not location_str:
            LOGGER.warning(f"Missing location for {circ_id}")
            return None

        location_data = helpers.parse_location_field(location_str)
        if not location_data:
            LOGGER.warning(f"Could not parse location for {circ_id}: {location_str}")
            return None

        # Get sequence from FASTA file
        sequence = sequences.get(circ_id, '')
        if not sequence:
            LOGGER.warning(f"No sequence found for {circ_id} in FASTA file")
            return None

        # Ensure sequence is DNA (T not U)
        sequence = sequence.upper().replace('U', 'T')

        # Get gene name from gene_Ensembl or gene_Refseq
        gene_name = row.get('gene_Refseq') or row.get('gene_Ensembl')
        if gene_name and not isinstance(gene_name, str):
            gene_name = str(gene_name)
        if gene_name:
            gene_name = gene_name.strip()

        # Parse exon positions from circname if available
        # Format: "circnhr-102(1-7)" indicates exons 1-7
        exons_list = []
        # Note: Exon positions would need to be derived from genome annotation
        # or additional data files - not directly available in this format

        # Build the Entry
        entry = Entry(
            primary_id=helpers.primary_id(circ_id),
            accession=helpers.accession(circ_id, location_data),
            ncbi_tax_id=taxid,
            database="CIRCPEDIA",
            sequence=sequence,
            regions=helpers.regions(
                location_data,
                location_data['strand'],
                exons_list,
                assembly_id
            ),
            rna_type=CIRCULAR_RNA_SO_TERM,
            url=helpers.url(circ_id),
            seq_version="1",
            note_data=helpers.note_data_from_tsv(row, location_data),
            chromosome=location_data['chromosome'],
            species=phy.species(taxid, taxonomy),
            common_name=phy.common_name(taxid, taxonomy),
            lineage=phy.lineage(taxid, taxonomy),
            gene=gene_name,
            product=helpers.product_from_gene(gene_name),
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
        LOGGER.error(f"Error parsing row for {circ_id if 'circ_id' in locals() else 'unknown'}: {e}", exc_info=True)
        return None


def parse(
    annotation_file: ty.Union[str, Path],
    fasta_file: ty.Union[str, Path],
    taxonomy_file: Path,
    assembly_id: ty.Optional[str] = None,
) -> ty.Iterable[Entry]:
    """
    Parse CIRCpedia V3 annotation and sequence files.

    Args:
        annotation_file: Path to CIRCpedia TSV annotation file
        fasta_file: Path to CIRCpedia FASTA sequence file
        taxonomy_file: Path to taxonomy database file
        assembly_id: Genome assembly ID (e.g., GRCh38, GRCm39)

    Yields:
        Entry objects for each circular RNA
    """
    LOGGER.info(f"Parsing CIRCpedia data from {annotation_file} and {fasta_file}")

    # Load taxonomy data
    taxonomy = SqliteDict(str(taxonomy_file))

    try:
        # Load sequences from FASTA
        sequences = load_fasta_sequences(fasta_file)

        # Read TSV annotation file with polars
        df = pl.read_csv(
            annotation_file,
            separator='\t',
            infer_schema_length=10000,
            ignore_errors=True,
        )

        LOGGER.info(f"Read {len(df)} rows from annotation file")
        LOGGER.info(f"Columns: {df.columns}")

        # Check for required columns
        if 'circID' not in df.columns and 'circid' not in df.columns:
            raise ValueError("Missing required column: circID")
        if 'Location' not in df.columns and 'location' not in df.columns:
            raise ValueError("Missing required column: Location")
        if 'species' not in df.columns:
            raise ValueError("Missing required column: species")

        # Process each row
        parsed_count = 0
        error_count = 0

        for row in df.iter_rows(named=True):
            entry = parse_tsv_row(row, sequences, taxonomy, assembly_id)
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
        LOGGER.error(f"Error parsing CIRCpedia files: {e}", exc_info=True)
        raise
    finally:
        taxonomy.close()
