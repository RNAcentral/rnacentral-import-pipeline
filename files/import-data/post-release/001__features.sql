\timing

BEGIN;

-- Drop indexes to speed up bulk insert
DROP INDEX IF EXISTS rnacen.ix_rnc_sequence_features__upi_taxid_name;
DROP INDEX IF EXISTS rnacen.ix_rnx_sequence_features_upi;

create index ix_load_rnc_sequence_features__accession on load_rnc_sequence_features(accession);

INSERT INTO rnc_sequence_features (
  upi,
  taxid,
  accession,
  start,
  stop,
  feature_name,
  metadata
) (
select distinct
  xref.upi,
  load.taxid,
  load.accession,
  load.start,
  load.stop,
  load.feature_name,
  load.metadata
from load_rnc_sequence_features load
join xref on xref.ac = load.accession
) ON CONFLICT (upi, taxid, accession, start, stop, feature_name) DO NOTHING;

drop table load_rnc_sequence_features;

-- Recreate indexes
CREATE INDEX ix_rnc_sequence_features__upi_taxid_name ON rnacen.rnc_sequence_features USING btree (upi, taxid, feature_name);
CREATE INDEX ix_rnx_sequence_features_upi ON rnacen.rnc_sequence_features USING btree (upi);

COMMIT;
