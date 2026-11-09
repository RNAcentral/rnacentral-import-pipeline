-- Sequence clustering (rnacentral-mmseqs-clustering): tier 1 absorbs fragments/variants into a
-- longer parent at 95% identity, tier 2 collapses near-duplicates at 90%. Keyed on URS, not
-- URS_taxid: the clustering is on sequence alone, so every taxon of a URS shares its flags.

-- Staging tables, filled by load-parquet from `rnac clustering parse`.
CREATE UNLOGGED TABLE IF NOT EXISTS load_sequence_cluster (
  urs text NOT NULL,
  cluster_95 text NOT NULL,
  cluster_90 text NOT NULL
);

CREATE UNLOGGED TABLE IF NOT EXISTS load_sequence_absorption (
  urs text NOT NULL,
  parent_urs text NOT NULL,
  absorbed_as text NOT NULL,
  identity real
);

-- One row per clustered sequence. Cluster ids are not stable between clustering runs, so this
-- is replaced wholesale on every load (see post-load.sql), never upserted.
CREATE TABLE IF NOT EXISTS rnc_sequence_cluster (
  urs text PRIMARY KEY REFERENCES rna(urs),
  -- Tier-1 parent. A sequence is its own parent unless something absorbed it.
  cluster_95 text NOT NULL REFERENCES rnc_sequence_cluster(urs) DEFERRABLE INITIALLY DEFERRED,
  -- Tier-2 representative, composed through the tier-1 parent. Always a tier-1 parent itself.
  cluster_90 text NOT NULL REFERENCES rnc_sequence_cluster(urs) DEFERRABLE INITIALLY DEFERRED,
  -- 'fragment': ratio to parent < 0.8; 'variant': near full length (ratio can exceed 1).
  absorbed_as text CHECK (absorbed_as IN ('fragment', 'variant')),
  -- Identity to the tier-1 parent; NULL when mmseqs align found no hit.
  identity_to_cluster_95 real,
  -- Never NULL for a clustered sequence; exports treat a missing row as "not clustered".
  is_fragment boolean GENERATED ALWAYS AS (absorbed_as IS NOT DISTINCT FROM 'fragment') STORED,
  in_deduplicated_set boolean GENERATED ALWAYS AS (urs = cluster_90) STORED,

  -- Absorbed exactly when it is not its own parent.
  CHECK ((urs = cluster_95) = (absorbed_as IS NULL)),
  CHECK (absorbed_as IS NOT NULL OR identity_to_cluster_95 IS NULL),
  -- Only tier-1 parents can be tier-2 representatives (the cross-row half is in post-load.sql).
  CHECK (urs <> cluster_90 OR urs = cluster_95)
);

CREATE INDEX IF NOT EXISTS ix__sequence_cluster__cluster_95 ON rnc_sequence_cluster(cluster_95);
CREATE INDEX IF NOT EXISTS ix__sequence_cluster__cluster_90 ON rnc_sequence_cluster(cluster_90);
-- The deduplicated set is a small fraction of the table and is what the exports read.
CREATE INDEX IF NOT EXISTS ix__sequence_cluster__dedup ON rnc_sequence_cluster(urs) WHERE in_deduplicated_set;
