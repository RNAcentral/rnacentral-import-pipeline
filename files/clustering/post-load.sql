-- Replace the whole clustering. Cluster ids change between runs, so merging row by row would
-- leave the old and new clusterings interleaved.
BEGIN;

-- The reps have to be present for the deferred self-references, and every absorbed member
-- must hang off the parent the clustering put it under; both mean a mismatched input pair.
DO $$
DECLARE
  bad bigint;
BEGIN
  SELECT count(*) INTO bad
  FROM load_sequence_absorption a
  LEFT JOIN load_sequence_cluster c ON c.urs = a.urs
  WHERE c.urs IS NULL OR c.cluster_95 <> a.parent_urs;
  IF bad > 0 THEN
    RAISE EXCEPTION '% absorbed members are missing from clusters.parquet or disagree with it on their tier-1 parent', bad;
  END IF;

  SELECT count(*) INTO bad
  FROM load_sequence_cluster c
  JOIN load_sequence_cluster p ON p.urs = c.cluster_90
  WHERE p.cluster_95 <> p.urs;
  IF bad > 0 THEN
    RAISE EXCEPTION '% sequences have a tier-2 representative that is not a tier-1 parent', bad;
  END IF;
END $$;

TRUNCATE rnc_sequence_cluster;

INSERT INTO rnc_sequence_cluster (
  urs,
  cluster_95,
  cluster_90,
  absorbed_as,
  identity_to_cluster_95
) (
SELECT
  c.urs,
  c.cluster_95,
  c.cluster_90,
  a.absorbed_as,
  a.identity
FROM load_sequence_cluster c
LEFT JOIN load_sequence_absorption a ON a.urs = c.urs
);

ANALYZE rnc_sequence_cluster;

TRUNCATE load_sequence_cluster, load_sequence_absorption;

COMMIT;
