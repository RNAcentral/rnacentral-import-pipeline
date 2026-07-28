-- One-off cleanup: drop redundant/unused indexes across the big rnacen tables.
-- Every entry is backed by pg_stat_user_indexes over a 14-day production window
-- (server up since 2026-06-29, stats never reset) - measured, not assumed.
--
-- Run AFTER the current load finishes. CONCURRENTLY = light lock, safe on a live
-- DB, must NOT run inside a transaction block:
--   psql "$PGDATABASE" -f files/drop-redundant-indexes.sql
--
-- Two of these (idx_rnc_sequence_regions_id, ix_rnc_rna_precomputed__upi_taxid_last_release)
-- are ALSO no longer recreated by the post-release pipeline (see
-- 000__populate_precompute.sql / 001__regions.sql); this drops the copies already
-- living in the DB. The rest are live-DB-only artifacts, so the drop is permanent.
--
--   table                  drop                                       scans/14d size   kept twin/cover
--   rnc_reference_map      ix_rnc_references__reference_id_accesion   0         54 GB  rnc_references_map$reference_id (1544)
--   rnc_reference_map      rnc_reference_map$reference_id             0         5.2GB  rnc_references_map$reference_id (dup)
--   rnc_reference_map      rnc_references_map$accession               391       27 GB  UNIQUE(accession,reference_id) (100M)
--   rnc_accession_active   idx__accession                            0         9.7GB  pkey (exact dup)
--   rnc_rna_precomputed    ix_..._upi_taxid_last_release             0         11 GB  upi_idx (upi,taxid, 423M)
--   rnc_sequence_regions   idx_rnc_sequence_regions_id               7784      1.8GB  pkey (exact dup, id)
--   rna                    idx_rna_upi                               0         1.7GB  rna_pkey (exact dup, 33.8M)
--   r2dt_results           ix_layout__urs                            79        1.2GB  un_layout__urs (UNIQUE)

-- Guard: refuse to run unless every fallback index we rely on still exists.
DO $$
DECLARE need text[] := ARRAY[
    'rnacen."rnc_references_map$accession$reference_id"',
    'rnacen."rnc_references_map$reference_id"',
    'rnacen.rnc_accession_active_pkey',
    'rnacen.rnc_rna_precomputed_upi_idx',
    'rnacen.rnc_sequence_regions_pkey',
    'rnacen.rna_pkey',
    'rnacen.un_layout__urs' ];
  r text;
BEGIN
  FOREACH r IN ARRAY need LOOP
    IF to_regclass(r) IS NULL THEN
      RAISE EXCEPTION 'fallback index % missing - aborting, dropping nothing', r;
    END IF;
  END LOOP;
  RAISE NOTICE 'all fallback indexes present - safe to drop redundant ones';
END $$;

DROP INDEX CONCURRENTLY IF EXISTS rnacen.ix_rnc_references__reference_id_accesion;         -- 54 GB, 0 scans
DROP INDEX CONCURRENTLY IF EXISTS rnacen."rnc_reference_map$reference_id";                 -- 5.2 GB, 0 scans, dup
DROP INDEX CONCURRENTLY IF EXISTS rnacen."rnc_references_map$accession";                   -- 27 GB, 391 scans
DROP INDEX CONCURRENTLY IF EXISTS rnacen.idx__accession;                                   -- 9.7 GB, 0 scans, dup of pkey
DROP INDEX CONCURRENTLY IF EXISTS rnacen.ix_rnc_rna_precomputed__upi_taxid_last_release;   -- 11 GB, 0 scans
DROP INDEX CONCURRENTLY IF EXISTS rnacen.idx_rnc_sequence_regions_id;                      -- 1.8 GB, dup of pkey
DROP INDEX CONCURRENTLY IF EXISTS rnacen.idx_rna_upi;                                      -- 1.7 GB, 0 scans, dup of pkey
DROP INDEX CONCURRENTLY IF EXISTS rnacen.ix_layout__urs;                                   -- 1.2 GB, dup of un_layout__urs

ANALYZE rnacen.rnc_reference_map;
ANALYZE rnacen.rnc_accession_active;
ANALYZE rnacen.rnc_rna_precomputed;
ANALYZE rnacen.rnc_sequence_regions;
ANALYZE rnacen.rna;
ANALYZE rnacen.r2dt_results;
