# -*- coding: utf-8 -*-

"""
Arrow schemas for the Parquet interchange format.

Each schema below is the canonical contract for one logical table produced by
``ParquetEntryWriter``. Column order matches the tuples yielded by the
corresponding ``Entry.write_*()`` method (or helper ``writeable()`` method on a
nested object). Reordering either side without updating the other will break
batch construction loudly — which is the whole point of declaring these.

All fields are currently ``nullable=True``. A follow-up pass will tighten NOT
NULL constraints once the canonical list is lifted out of Postgres.
"""

from __future__ import annotations

import pyarrow as pa

# ---------------------------------------------------------------------------
# Entry.write_ac_info() -> accessions.parquet
# Source of truth for column names: files/import-data/load/accessions.ctl
# NOT NULLs reflect destination rnc_accessions constraints (only 'accession').
# Sequence-related fields (short/long_sequences) are held for the
# update_rnc_accessions() stored-procedure audit.
ACCESSIONS = pa.schema(
    [
        pa.field("accession", pa.string(), nullable=False),
        pa.field("parent_ac", pa.string()),
        pa.field("seq_version", pa.string()),
        pa.field("feature_start", pa.int64()),
        pa.field("feature_end", pa.int64()),
        pa.field("feature_name", pa.string()),
        pa.field("non_coding_id", pa.string()),
        pa.field("database", pa.string()),
        pa.field("external_id", pa.string()),
        pa.field("optional_id", pa.string()),
        pa.field("description", pa.string()),
        pa.field("organelle", pa.string()),
        pa.field("chromosome", pa.string()),
        pa.field("function", pa.string()),
        pa.field("gene", pa.string()),
        pa.field("gene_synonym", pa.string()),
        pa.field("inference", pa.string()),
        pa.field("locus_tag", pa.string()),
        pa.field("mol_type", pa.string()),
        pa.field("ncrna_class", pa.string()),
        pa.field("note", pa.string()),
        pa.field("product", pa.string()),
        pa.field("standard_name", pa.string()),
        pa.field("db_xref", pa.string()),
        pa.field("so_term", pa.string()),
        pa.field("url", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_seq_short() / write_seq_long() -> short/long_sequences.parquet
# Source of truth: files/import-data/load/short-sequences.ctl, long-sequences.ctl
# Same 9-column shape; split purely by sequence length (<=4000 vs >4000) for
# legacy Postgres schema reasons (ex-Oracle VARCHAR2 limit).
# TODO: tighten NOT NULLs (all columns in practice)
# NOT NULLs per the downstream functions: they dedupe/lookup sequences by
# combinations of (database, ac, version, taxid, crc64, md5) and the sequence
# split by database id, so every one of those must be present. optional_id is
# the one genuinely optional field.
_SEQUENCES_FIELDS = [
    pa.field("crc64", pa.string(), nullable=False),
    pa.field("len", pa.int64(), nullable=False),
    pa.field("database", pa.string(), nullable=False),
    pa.field("ac", pa.string(), nullable=False),
    pa.field("optional_id", pa.string()),
    pa.field("version", pa.string(), nullable=False),
    pa.field("taxid", pa.int64(), nullable=False),
    pa.field("md5", pa.string(), nullable=False),
]

SHORT_SEQUENCES = pa.schema(
    [
        _SEQUENCES_FIELDS[0],
        _SEQUENCES_FIELDS[1],
        pa.field("seq_short", pa.string(), nullable=False),
    ]
    + _SEQUENCES_FIELDS[2:]
)

LONG_SEQUENCES = pa.schema(
    [
        _SEQUENCES_FIELDS[0],
        _SEQUENCES_FIELDS[1],
        pa.field("seq_long", pa.string(), nullable=False),
    ]
    + _SEQUENCES_FIELDS[2:]
)


# ---------------------------------------------------------------------------
# Entry.write_refs() -> references.parquet
# Source of truth: files/import-data/load/references.ctl
# Reference.writeable() yields: [md5, accession, authors, location, title, pmid, doi]
# NOT NULLs:
#   - md5: enforced on rnc_references (PK).
#   - accession: feeds rnc_reference_map.accession in the second INSERT of
#     rnc_update.update_literature_references(), which requires non-null.
REFERENCES = pa.schema(
    [
        pa.field("md5", pa.string(), nullable=False),
        pa.field("accession", pa.string(), nullable=False),
        pa.field("authors", pa.string()),
        pa.field("location", pa.string()),
        pa.field("title", pa.string()),
        pa.field("pmid", pa.int64()),
        pa.field("doi", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_ref_ids() -> ref_ids.parquet
# IdReference.writeable(accession) yields: [normalized_id, accession]
# No .ctl currently; kept for parity with the CSV writer's output.
REF_IDS = pa.schema(
    [
        pa.field("normalized_id", pa.string()),
        pa.field("accession", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_sequence_regions() -> regions.parquet
# Source: files/import-data/load/regions.ctl
# SequenceRegion.writeable(accession) yields one row per exon.
# NOT NULLs reflect rnc_sequence_regions (via 001__regions.sql):
#   - urs_taxid is derived by joining on accession, so accession must be non-null.
#   - region_name, chromosome, strand, exon_count are passed through directly.
#   - region_start/region_stop on the destination are min/max of exon_start/stop,
#     so the exon coordinates themselves must be non-null.
#   - assembly_id stays nullable (destination allows NULL).
REGIONS = pa.schema(
    [
        pa.field("accession", pa.string(), nullable=False),
        pa.field("region_name", pa.string(), nullable=False),
        pa.field("chromosome", pa.string(), nullable=False),
        pa.field("strand", pa.int32(), nullable=False),
        pa.field("assembly_id", pa.string()),
        pa.field("exon_count", pa.int32(), nullable=False),
        pa.field("exon_start", pa.int64(), nullable=False),
        pa.field("exon_stop", pa.int64(), nullable=False),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_secondary_structure() -> secondary_structure.parquet
# Source: files/import-data/load/secondary-structure.ctl
# SecondaryStructure.writeable(accession) yields: [accession, dot_bracket, md5]
SECONDARY_STRUCTURE = pa.schema(
    [
        pa.field("rnc_accession_id", pa.string()),
        pa.field("secondary_structure", pa.string()),
        pa.field("md5", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_related_sequences() -> related_sequences.parquet
# Source: files/import-data/load/related-sequences.ctl
# RelatedSequence.writeable(accession) yields:
#   [source_accession, target_accession, relationship_type, methods]
# `methods` is currently a Postgres-array literal string like '{"a","b"}'. Kept
# as pa.string() so behaviour matches the CSV path exactly; tighten to
# pa.list_(pa.string()) in a later pass when we know DuckDB's COPY is happy.
# NOT NULLs on rnc_related_sequences: source_accession, target_accession,
# relationship_type (source_urs_taxid is derived from source_accession via a
# join). methods stays nullable.
RELATED_SEQUENCES = pa.schema(
    [
        pa.field("source_accession", pa.string(), nullable=False),
        pa.field("target_accession", pa.string(), nullable=False),
        pa.field("relationship_type", pa.string(), nullable=False),
        pa.field("methods", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_sequence_features() -> features.parquet
# Source: files/import-data/load/features.ctl
# SequenceFeature.writeable(accession, taxid) yields:
#   [accession, taxid, start, stop, feature_name, metadata(json string)]
# NOT NULLs on rnc_sequence_features: upi (from accession join), start, stop,
# feature_name. taxid is NOT NULL on load_rnc_sequence_features itself.
# metadata is jsonb-typed and nullable.
FEATURES = pa.schema(
    [
        pa.field("accession", pa.string(), nullable=False),
        pa.field("taxid", pa.int64(), nullable=False),
        pa.field("start", pa.int64(), nullable=False),
        pa.field("stop", pa.int64(), nullable=False),
        pa.field("feature_name", pa.string(), nullable=False),
        pa.field("metadata", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_interactions() -> interactions.parquet
# Source: files/import-data/load/interactions.ctl
# Interaction.writeable() yields: [intact_id, urs_taxid, interacting_id, names(json string), taxid]
# Every column NOT NULL on both load_interactions and rnc_interactions.
INTERACTIONS = pa.schema(
    [
        pa.field("intact_id", pa.string(), nullable=False),
        pa.field("urs_taxid", pa.string(), nullable=False),
        pa.field("interacting_id", pa.string(), nullable=False),
        pa.field("names", pa.string(), nullable=False),
        pa.field("taxid", pa.int64(), nullable=False),
    ]
)


# ---------------------------------------------------------------------------
# Entry.write_ontology_terms() -> terms.parquet
# 1-column listing of SO term IDs referenced by entries + evidence codes from
# GO annotations (see writers.py). Downstream join reconstructs full ontology
# metadata; this file is the set of terms we referenced, not the ontology.
TERMS = pa.schema(
    [
        pa.field("ontology_term_id", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# GoTermAnnotation.writeable() -> go_annotations.parquet
# Source: files/import-data/load/go-annotations.ctl
# NOT NULLs on go_term_annotations: rna_id, qualifier, ontology_term_id,
# evidence_code. assigned_by and extensions stay nullable (destination allows
# NULL even though un_go_term_annotations uniqueness includes assigned_by).
GO_ANNOTATIONS = pa.schema(
    [
        pa.field("rna_id", pa.string(), nullable=False),
        pa.field("qualifier", pa.string(), nullable=False),
        pa.field("ontology_term_id", pa.string(), nullable=False),
        pa.field("evidence_code", pa.string(), nullable=False),
        pa.field("extensions", pa.string()),
        pa.field("assigned_by", pa.string()),
    ]
)


# ---------------------------------------------------------------------------
# GoTermAnnotation.writeable_publication_mappings() -> go_publication_mappings.parquet
# Source: files/import-data/load/go-publication-mappings.ctl
# go_term_publication_map has two NOT NULL columns (go_term_annotation_id,
# reference_id), both resolved by joins: the annotation FK joins on
# (rna_id, qualifier, ontology_term_id, evidence_code, assigned_by) against
# go_term_annotations, and reference_id joins on pubmed_id. All six parser
# fields therefore must be non-null for the post-release INSERT to succeed.
GO_PUBLICATION_MAPPINGS = pa.schema(
    [
        pa.field("rna_id", pa.string(), nullable=False),
        pa.field("qualifier", pa.string(), nullable=False),
        pa.field("ontology_term_id", pa.string(), nullable=False),
        pa.field("assigned_by", pa.string(), nullable=False),
        pa.field("evidence_code", pa.string(), nullable=False),
        pa.field("pubmed_id", pa.string(), nullable=False),
    ]
)


# ---------------------------------------------------------------------------
# rfam.infernal_results.as_csv() -> hits.parquet (rfam-scan workflow)
# Source: files/rfam-scan/load.ctl + files/schema/create_load.sql.
# Column names match the destination ``load_rfam_model_hits`` table directly
# (parquet col names must match destination col names; bin/load-parquet has
# no equivalent of pgloader's TARGET COLUMNS re-mapping). The ``strand``
# column the writer historically emitted is dropped here — pgloader was
# already discarding it via the ctl's narrower TARGET COLUMNS list.
# All nine columns are NOT NULL on the destination table.
RFAM_HITS = pa.schema(
    [
        pa.field("upi", pa.string(), nullable=False),
        pa.field("sequence_start", pa.int64(), nullable=False),
        pa.field("sequence_stop", pa.int64(), nullable=False),
        pa.field("rfam_model_id", pa.string(), nullable=False),
        pa.field("model_start", pa.int64(), nullable=False),
        pa.field("model_stop", pa.int64(), nullable=False),
        pa.field("overlap", pa.string(), nullable=False),
        pa.field("e_value", pa.float64(), nullable=False),
        pa.field("score", pa.float64(), nullable=False),
    ]
)


# ---------------------------------------------------------------------------
# Logical name -> Postgres staging table name. Lifted from the INTO clauses
# of files/import-data/load/*.ctl so the Parquet load path mirrors pgloader's
# table targets exactly. Used by bin/load-parquet when it receives a logical
# name (matching a key in ENTRY_WRITER_SCHEMAS) instead of an explicit table.
#
# Note: ref_ids has no pgloader ctl today; it's parser output the CSV path
# never loaded. Mapping to None signals "skip load" rather than a typo.
# terms.parquet is a 1-column listing used downstream for joins and doesn't
# land in load_ontology_terms (which has a different 4-column shape) either.
ENTRY_WRITER_LOAD_TABLES: "dict[str, str | None]" = {
    "accessions": "load_rnc_accessions",
    "short_sequences": "load_rnacentral_all",
    "long_sequences": "load_rnacentral_all",
    "references": "load_rnc_references",
    "ref_ids": None,
    "regions": "load_rnc_sequence_regions",
    "secondary_structure": "load_rnc_secondary_structure",
    "related_sequences": "load_rnc_related_sequences",
    "features": "load_rnc_sequence_features",
    "interactions": "load_interactions",
    "terms": None,
    "go_annotations": "load_go_term_annotations",
    "go_publication_mappings": "load_go_term_publication_map",
}


# ---------------------------------------------------------------------------
# Mapping used by ParquetEntryWriter. Keys match EntryWriter attr names so the
# two writers stay structurally aligned.
ENTRY_WRITER_SCHEMAS = {
    "accessions": ACCESSIONS,
    "short_sequences": SHORT_SEQUENCES,
    "long_sequences": LONG_SEQUENCES,
    "references": REFERENCES,
    "ref_ids": REF_IDS,
    "regions": REGIONS,
    "secondary_structure": SECONDARY_STRUCTURE,
    "related_sequences": RELATED_SEQUENCES,
    "features": FEATURES,
    "interactions": INTERACTIONS,
    "terms": TERMS,
    "go_annotations": GO_ANNOTATIONS,
    "go_publication_mappings": GO_PUBLICATION_MAPPINGS,
}
