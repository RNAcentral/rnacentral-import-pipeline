-- Merge fresh editing-event features from load_rediportal_features into
-- rnc_sequence_features. Lifted from the AFTER LOAD DO body of
-- files/rediportal/load.ctl. Run by bin/load-parquet --post-load after the
-- staging table is populated; the entire script executes in a single
-- transaction so a failure leaves rnc_sequence_features untouched.
--
-- IMPORTANT: the DELETE is unconditional on feature_name, NOT joined to the
-- staging table. An empty load would therefore wipe all rna_editing_event
-- rows from production. ``bin/load-parquet`` skips the post-load step when
-- 0 rows were loaded (unless --allow-empty is passed) — that guard is
-- essential here.

DELETE FROM rnc_sequence_features WHERE feature_name = 'rna_editing_event';

INSERT INTO rnc_sequence_features (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata,
    feature_provider
)
SELECT DISTINCT
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata,
    feature_provider
FROM load_rediportal_features;
