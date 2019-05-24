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
) ON CONFLICT (un_rnc_sequence_features__upi_taxid_accession_start_stop_feature_name) DO NOTHING;

drop table load_rnc_sequence_features;
