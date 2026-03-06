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

import logging
import re
import typing as ty
from urllib.parse import quote

from rnacentral_pipeline.databases.data import (
    CoordinateSystem,
    Exon,
    Reference,
    SequenceRegion,
)
from rnacentral_pipeline.databases.data.features import SequenceFeature
from rnacentral_pipeline.databases.data.related import RelatedEvidence, RelatedSequence
from rnacentral_pipeline.databases.helpers import phylogeny as phy

LOGGER = logging.getLogger(__name__)

# CIRCpedia base URL
CIRCPEDIA_BASE_URL = "https://bits.fudan.edu.cn/circpediav3"


def get_ncbi_taxid(species_str: str) -> ty.Optional[int]:
    """
    Get NCBI taxonomy ID from species name.

    Args:
        species_str: Species string

    Returns:
        NCBI taxonomy ID or None if not found
    """
    try:
        # Try to get taxid from scientific name
        taxid = phy.taxid(species_str)
        if taxid:
            return taxid
    except (phy.UnknownTaxonId, phy.FailedTaxonId, Exception):
        raise ValueError(f"Failed to lookup taxid for {species_str}")

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


def parse_location_field(location_str: str) -> ty.Optional[ty.Dict[str, ty.Any]]:
    """
    Parse CIRCpedia V3 Location field which combines chromosome, coordinates, and strand.

    Expected format: "V:15874634-15876408(-)" or "chr1:100-200(+)"
    Components:
        - Chromosome: V, chr1, I, II, X, etc.
        - Start:End coordinates (1-based, fully-closed)
        - Strand: (+) or (-)

    Args:
        location_str: Location string from CIRCpedia V3

    Returns:
        Dict with chromosome, start, end, strand or None if parsing fails
    """
    if not location_str:
        return None

    # Match pattern: chromosome:start-end(strand)
    # Examples: "V:15874634-15876408(-)", "chr1:100-200(+)", "III:1000-2000(+)"
    match = re.match(
        r"^(?:chr)?([^:]+):(\d+)-(\d+)\(([+-])\)$", str(location_str).strip()
    )
    if not match:
        LOGGER.warning(f"Could not parse location field: {location_str}")
        return None

    chromosome, start, end, strand_char = match.groups()
    # Convert strand character to integer
    strand = 1 if strand_char == "+" else -1

    return {
        "chromosome": chromosome,
        "start": int(start),
        "end": int(end),
        "strand": strand,
    }


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

    Uses CIRCpedia's own identifier prefixed with database name, with location
    to ensure uniqueness when a circRNA appears in multiple locations.

    IMPORTANT: This is NOT the URS identifier. RNAcentral automatically
    generates URS identifiers from sequence hashes. This accession is
    the external database's own identifier that links back to CIRCpedia.

    Args:
        circ_id: CIRCpedia circRNA ID
        location: Location dict with chromosome, start, end

    Returns:
        Accession string prefixed with database name (e.g., "CIRCPEDIA:hsa_circ_0001_1:100-200")
    """
    # Prefix with database name to avoid ID conflicts across databases
    # Include location to ensure uniqueness when a circRNA appears at multiple locations
    location_str = f"{location['chromosome']}:{location['start']}-{location['end']}"
    return f"CIRCPEDIA:{circ_id}|{location_str}"


def url(circ_id: str) -> str:
    """
    Generate URL for circRNA in CIRCpedia V3.

    Returns a direct link to the circRNA detail page.
    Format: https://bits.fudan.edu.cn/circpediav3/circrna/{circ_id}

    Args:
        circ_id: CIRCpedia circRNA ID

    Returns:
        URL string to the circRNA detail page
    """
    # Encode the ID for URL safety
    encoded_id = quote(str(circ_id))
    return f"{CIRCPEDIA_BASE_URL}/circrna/{encoded_id}"


def regions(
    location: ty.Dict[str, ty.Any],
    assembly_id: ty.Optional[str] = None,
) -> ty.List[SequenceRegion]:
    """
    Build sequence regions from location data.

    CIRCpedia uses standard genomic coordinates (chr:start-end format)
    which follow the 1-based, fully-closed coordinate system (same as
    GFF/GTF format). Both start and end positions are inclusive.

    Args:
        location: Location dict with chromosome, start, end, strand
        assembly_id: Assembly ID (e.g., GRCh38)

    Returns:
        List of SequenceRegion objects
    """
    # Create a single exon spanning the entire circular RNA region
    # Note: CIRCpedia doesn't provide individual exon coordinates
    exons = [Exon(start=location["start"], stop=location["end"])]

    return [
        SequenceRegion(
            chromosome=location["chromosome"],
            strand=location["strand"],
            exons=exons,
            assembly_id=assembly_id or "unknown",
            coordinate_system=CoordinateSystem.one_based(),
        )
    ]


def note_data_from_tsv(
    row: ty.Dict[str, ty.Any], location: ty.Dict[str, ty.Any]
) -> ty.Dict[str, ty.Any]:
    """
    Build note data dictionary from CIRCpedia V3 TSV row.

    Extracts metadata fields from the actual TSV columns:
    circID, Location, circname, editing_site, subcell_location, length,
    DIS3_signal, DIS3_motif, Orthology, TGS, species, gene_Ensembl,
    gene_Refseq, transcript_Ensembl, transcript_Refseq

    Args:
        row: Dictionary representing TSV row
        location: Location dict with chromosome, start, end

    Returns:
        Note data dictionary
    """
    notes = {}

    # Add genomic location
    notes[
        "genomic_location"
    ] = f"{location['chromosome']}:{location['start']}-{location['end']}"

    # Add circular RNA name if available
    if row.get("circname"):
        notes["circname"] = str(row["circname"])

    # Add length
    if row.get("length"):
        try:
            notes["length"] = int(row["length"])
        except (ValueError, TypeError):
            pass

    # Add subcellular localization
    if row.get("subcell_location"):
        notes["subcell_location"] = str(row["subcell_location"])

    # Add editing sites if available (filter 'none')
    if row.get("editing_site") and set(row["editing_site"].split(";")) != {"none"}:
        notes["editing_site"] = str(row["editing_site"])

    # Add DIS3 degradation signals (filter 'none')
    if row.get("DIS3_signal") and row["DIS3_signal"] != "none":
        notes["DIS3_signal"] = str(row["DIS3_signal"])

    # Add DIS3 motif if available (filter 'none')
    if row.get("DIS3_motif") and row["DIS3_motif"] != "none":
        notes["DIS3_motif"] = str(row["DIS3_motif"])

    # Add orthology information (filter 'none')
    if row.get("Orthology") and row["Orthology"] != "none":
        notes["orthology"] = str(row["Orthology"])

    # Add TGS (third-generation sequencing) support (filter 'none')
    if row.get("TGS") and row["TGS"] != "none":
        notes["TGS_support"] = str(row["TGS"])

    # Add transcript IDs (filter 'NA')
    if row.get("transcript_Ensembl") and row["transcript_Ensembl"] != "NA":
        notes["transcript_ensembl"] = str(row["transcript_Ensembl"])
    if row.get("transcript_Refseq") and row["transcript_Refseq"] != "NA":
        notes["transcript_refseq"] = str(row["transcript_Refseq"])

    return notes


def product_from_gene(gene_name: ty.Optional[str]) -> str:
    """
    Generate product description from gene name.

    Args:
        gene_name: Gene name

    Returns:
        Product description
    """
    if gene_name:
        return f"{gene_name} circular RNA"
    return "circular RNA"


def description(
    taxid: int,
    gene_name: ty.Optional[str],
) -> str:
    """
    Generate description for the entry.

    Args:
        taxid: NCBI taxonomy ID
        gene_name: Gene name

    Returns:
        Description string
    """
    try:
        species_name = phy.species(taxid)
        common = phy.common_name(taxid)

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


def dis3_features(dis3_motif: str) -> ty.List[SequenceFeature]:
    """
    Parse DIS3 degradation signal motifs into SequenceFeature objects.

    Args:
        dis3_motif: DIS3_motif column value, e.g.
            "273-279:GTCCTGC:M2-GCTCTGC|308-314:GTTGCTG:optiM-TTTGCTG"
            Multiple motif matches at one signal are comma-separated within a
            pipe-delimited group.

    Returns:
        List of SequenceFeature objects (empty if value is 'none' or empty)
    """
    if not dis3_motif or dis3_motif.strip().lower() == "none":
        return []

    features = []
    for group in dis3_motif.split("|"):
        for motif_entry in group.split(","):
            motif_entry = motif_entry.strip()
            if not motif_entry:
                continue

            # Format: start-end:observed_seq:class-consensus
            parts = motif_entry.split(":")
            if len(parts) != 3:
                LOGGER.warning(f"Unexpected DIS3 motif format: {motif_entry}")
                continue

            coord_str, observed_seq, class_consensus = parts

            # Parse coordinates (1-based closed → 0-based half-open)
            try:
                start_str, end_str = coord_str.split("-")
                start = int(start_str) - 1  # convert to 0-based
                end = int(end_str)  # half-open
            except (ValueError, TypeError):
                LOGGER.warning(f"Could not parse DIS3 coordinates: {coord_str}")
                continue

            # Split class-consensus: first token is class, rest is consensus
            # e.g. "M2-GCTCTGC" → class="M2", consensus="GCTCTGC"
            # e.g. "optiM-TTTGCTG" → class="optiM", consensus="TTTGCTG"
            cc_parts = class_consensus.split("-", 1)
            motif_class = cc_parts[0] if len(cc_parts) >= 1 else class_consensus
            consensus = cc_parts[1] if len(cc_parts) >= 2 else ""

            features.append(
                SequenceFeature(
                    name="dis3_degradation_signal",
                    feature_type="dis3_degradation_signal",
                    location=[start, end],
                    sequence=observed_seq.upper(),
                    provider="CIRCPEDIA",
                    metadata={"motif_class": motif_class, "consensus": consensus},
                )
            )

    return features


def ortholog_sequences(orthology: str) -> ty.List[RelatedSequence]:
    """
    Parse ortholog circRNA IDs from the Orthology column into RelatedSequence objects.

    Args:
        orthology: Comma-separated ortholog circRNA IDs, or "none" if no orthologs

    Returns:
        List of RelatedSequence objects
    """
    if not orthology or orthology.strip().lower() == "none":
        return []

    related = []
    for ortholog_id in orthology.split(","):
        ortholog_id = ortholog_id.strip()
        if not ortholog_id:
            continue
        related.append(
            RelatedSequence(
                sequence_id=f"CIRCPEDIA:{ortholog_id}",
                relationship="homolog",
                evidence=RelatedEvidence(methods=["liftover - provided"]),
            )
        )
    return related


def references() -> ty.List[Reference]:
    """
    Get references for CIRCpedia V3.

    Returns:
        List of Reference objects
    """
    return [
        Reference(
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
