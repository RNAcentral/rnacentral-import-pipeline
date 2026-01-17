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
Helper functions for parsing CIRCpedia V3 data.
"""

import hashlib
import logging
import re
import typing as ty
from urllib.parse import quote

import polars as pl
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.data import (
    Entry,
    Exon,
    IdReference,
    Reference,
    SequenceRegion,
)
from rnacentral_pipeline.databases.helpers import phylogeny as phy

LOGGER = logging.getLogger(__name__)

# CIRCpedia base URL
CIRCPEDIA_BASE_URL = "https://bits.fudan.edu.cn/circpediav3"


def parse_species_name(species_str: str) -> ty.Tuple[str, str]:
    """
    Parse species name from CIRCpedia format.

    Args:
        species_str: Species string (e.g., "Homo sapiens", "Mus musculus")

    Returns:
        Tuple of (genus, species)
    """
    parts = species_str.strip().split()
    if len(parts) >= 2:
        return parts[0], parts[1]
    return species_str, ""


def get_ncbi_taxid(taxonomy: SqliteDict, species_str: str) -> ty.Optional[int]:
    """
    Get NCBI taxonomy ID from species name.

    Args:
        taxonomy: SqliteDict with taxonomy data
        species_str: Species string

    Returns:
        NCBI taxonomy ID or None if not found
    """
    try:
        # Try to get taxid from scientific name
        taxid = phy.name_to_taxid(species_str, taxonomy)
        if taxid:
            return taxid
    except (phy.UnknownTaxonId, phy.FailedTaxonId, Exception):
        pass

    # Common species mapping as fallback
    species_map = {
        "homo sapiens": 9606,
        "mus musculus": 10090,
        "rattus norvegicus": 10116,
        "danio rerio": 7955,
        "drosophila melanogaster": 7227,
        "caenorhabditis elegans": 6239,
        "saccharomyces cerevisiae": 4932,
        "arabidopsis thaliana": 3702,
        "gallus gallus": 9031,
        "pan troglodytes": 9598,
        "macaca mulatta": 9544,
        "canis familiaris": 9615,
        "bos taurus": 9913,
        "sus scrofa": 9823,
        "xenopus tropicalis": 8364,
        "oryzias latipes": 8090,
        "takifugu rubripes": 31033,
        "strongylocentrotus purpuratus": 7668,
        "ciona intestinalis": 7719,
        "anopheles gambiae": 7165,
    }

    species_lower = species_str.lower().strip()
    if species_lower in species_map:
        return species_map[species_lower]

    LOGGER.warning(f"Could not find taxid for species: {species_str}")
    return None


def parse_genomic_location(location_str: str) -> ty.Optional[ty.Dict[str, ty.Any]]:
    """
    Parse genomic location string from CIRCpedia.

    Expected formats:
        - chr1:12345-67890
        - 1:12345-67890
        - chrX:12345-67890

    Args:
        location_str: Location string

    Returns:
        Dict with chromosome, start, end or None if parsing fails
    """
    if not location_str or pd.isna(location_str):
        return None

    # Remove 'chr' prefix if present for normalized chromosome name
    match = re.match(r'^(?:chr)?([^:]+):(\d+)-(\d+)$', str(location_str).strip())
    if not match:
        LOGGER.warning(f"Could not parse location: {location_str}")
        return None

    chromosome, start, end = match.groups()
    return {
        "chromosome": chromosome,
        "start": int(start),
        "end": int(end),
    }


def parse_exon_positions(exon_str: str) -> ty.List[ty.Tuple[int, int]]:
    """
    Parse exon positions from CIRCpedia format.

    Expected format: "start1-end1,start2-end2,..."

    Args:
        exon_str: Exon positions string

    Returns:
        List of (start, end) tuples
    """
    if not exon_str or pd.isna(exon_str):
        return []

    exons = []
    for exon_pos in str(exon_str).split(','):
        exon_pos = exon_pos.strip()
        if '-' in exon_pos:
            try:
                start, end = exon_pos.split('-')
                exons.append((int(start), int(end)))
            except (ValueError, AttributeError):
                LOGGER.warning(f"Could not parse exon position: {exon_pos}")

    return exons


def parse_strand(strand_str: str) -> int:
    """
    Parse strand information.

    Args:
        strand_str: Strand string ('+' or '-')

    Returns:
        1 for forward, -1 for reverse, 0 for unknown
    """
    if not strand_str or pd.isna(strand_str):
        return 0

    strand_str = str(strand_str).strip()
    if strand_str == '+':
        return 1
    elif strand_str == '-':
        return -1
    return 0


def primary_id(circ_id: str) -> str:
    """
    Generate primary ID for RNAcentral.

    Args:
        circ_id: CIRCpedia circRNA ID

    Returns:
        Primary ID string
    """
    return str(circ_id)


def accession(circ_id: str, location: ty.Dict[str, ty.Any]) -> str:
    """
    Generate accession for RNAcentral.

    This creates a unique accession based on circRNA ID and location.

    Args:
        circ_id: CIRCpedia circRNA ID
        location: Location dict with chromosome, start, end

    Returns:
        Accession string
    """
    # Create a unique accession based on ID and location
    location_str = f"{location['chromosome']}:{location['start']}-{location['end']}"
    combined = f"{circ_id}_{location_str}"

    # Use hash to create a unique but consistent accession
    hash_val = hashlib.md5(combined.encode()).hexdigest()[:16]
    return f"CIRCPEDIA_{hash_val}".upper()


def url(circ_id: str) -> str:
    """
    Generate URL for circRNA in CIRCpedia.

    Args:
        circ_id: CIRCpedia circRNA ID

    Returns:
        URL string
    """
    # Encode the ID for URL safety
    encoded_id = quote(str(circ_id))
    return f"{CIRCPEDIA_BASE_URL}/search?query={encoded_id}"


def seq_version(row: pl.Series) -> str:
    """
    Get sequence version (default to 1).

    Args:
        row: Polars row

    Returns:
        Version string
    """
    return "1"


def note_data(row: pl.Series, location: ty.Dict[str, ty.Any]) -> ty.Dict[str, ty.Any]:
    """
    Build note data dictionary with additional information.

    Args:
        row: Polars row with circRNA data
        location: Location dict

    Returns:
        Note data dictionary
    """
    notes = {}

    # Add expression data if available (FPM = Fragments Per Million)
    if "fpm" in row and row["fpm"] is not None and not pd.isna(row["fpm"]):
        notes["expression_fpm"] = float(row["fpm"])

    # Add cell line/tissue information
    if "cell_line" in row and row["cell_line"] is not None:
        notes["cell_line"] = str(row["cell_line"])

    if "tissue" in row and row["tissue"] is not None:
        notes["tissue"] = str(row["tissue"])

    # Add sequencing type
    if "sequencing_type" in row and row["sequencing_type"] is not None:
        notes["sequencing_type"] = str(row["sequencing_type"])

    # Add conservation information
    if "conservation" in row and row["conservation"] is not None:
        notes["conservation"] = str(row["conservation"])

    # Add RNase R enrichment
    if "rnase_r_enrichment" in row and row["rnase_r_enrichment"] is not None:
        notes["rnase_r_enrichment"] = float(row["rnase_r_enrichment"])

    # Add genomic location
    notes["genomic_location"] = f"{location['chromosome']}:{location['start']}-{location['end']}"

    return notes


def chromosome(location: ty.Dict[str, ty.Any]) -> str:
    """
    Get chromosome name.

    Args:
        location: Location dict

    Returns:
        Chromosome name
    """
    return location["chromosome"]


def regions(
    location: ty.Dict[str, ty.Any],
    strand: int,
    exons_list: ty.List[ty.Tuple[int, int]],
    assembly_id: ty.Optional[str] = None,
) -> ty.List[SequenceRegion]:
    """
    Build sequence regions from location and exon data.

    Args:
        location: Location dict with chromosome, start, end
        strand: Strand (1, -1, or 0)
        exons_list: List of (start, end) tuples for exons
        assembly_id: Assembly ID (e.g., GRCh38)

    Returns:
        List of SequenceRegion objects
    """
    # If no exons provided, create a single exon spanning the entire region
    if not exons_list:
        exons_list = [(location["start"], location["end"])]

    exons = [Exon(start=start, stop=end) for start, end in exons_list]

    return [
        SequenceRegion(
            chromosome=location["chromosome"],
            strand=strand,
            exons=exons,
            assembly_id=assembly_id or "unknown",
            coordinate_system=None,
        )
    ]


def gene(row: pl.Series) -> ty.Optional[str]:
    """
    Get gene name.

    Args:
        row: Polars row

    Returns:
        Gene name or None
    """
    if "gene" in row and row["gene"] is not None and not pd.isna(row["gene"]):
        return str(row["gene"])
    return None


def product(row: pl.Series, gene_name: ty.Optional[str]) -> str:
    """
    Generate product description.

    Args:
        row: Polars row
        gene_name: Gene name

    Returns:
        Product description
    """
    if gene_name:
        return f"{gene_name} circular RNA"
    return "circular RNA"


def description(
    taxonomy: SqliteDict,
    taxid: int,
    gene_name: ty.Optional[str],
) -> str:
    """
    Generate description for the entry.

    Args:
        taxonomy: SqliteDict with taxonomy data
        taxid: NCBI taxonomy ID
        gene_name: Gene name

    Returns:
        Description string
    """
    try:
        species_name = phy.species(taxid, taxonomy)
        common = phy.common_name(taxid, taxonomy)

        if common:
            prefix = f"{species_name} ({common})"
        else:
            prefix = species_name

        if gene_name:
            return f"{prefix} {gene_name} circular RNA"
        return f"{prefix} circular RNA"
    except (phy.UnknownTaxonId, phy.FailedTaxonId, Exception):
        if gene_name:
            return f"{gene_name} circular RNA"
        return "circular RNA"


def references() -> ty.List[Reference]:
    """
    Get references for CIRCpedia V3.

    Returns:
        List of Reference objects
    """
    return [
        IdReference(
            authors="Zhai SN, Zhang YY, Chen MH, Fu ZC, Chen LL, Ma XK, Yang L",
            location=(
                "CIRCpedia v3: an interactive database for circular RNA "
                "characterization and functional exploration"
            ),
            title=(
                "CIRCpedia v3: an interactive database for circular RNA "
                "characterization and functional exploration"
            ),
            pmid=None,  # Update when PMID is available
            doi="10.1093/nar/gkaf1039",
        )
    ]


# Import pandas for compatibility check
try:
    import pandas as pd
except ImportError:
    # Create a minimal mock for isna if pandas not available
    class pd:
        @staticmethod
        def isna(value):
            return value is None or (isinstance(value, float) and value != value)
