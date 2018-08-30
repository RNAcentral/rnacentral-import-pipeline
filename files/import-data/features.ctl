LOAD CSV
FROM ALL FILENAMES MATCHING ~<features*.csv>
HAVING FIELDS (
    accession,
    taxid,
    start,
    stop,
    feature_name,
    metadata
)
INTO {{PGDATABASE}}?load_rnc_sequence_features
TARGET COLUMNS (
    accession,
    taxid,
    start,
    stop,
    feature_name,
    metadata
)

WITH truncate,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_rnc_sequence_features;
$$,
$$
create table load_rnc_sequence_features (
    accession varchar(100) NOT NULL,
    taxid int not null,
    start int not null,
    stop int not null,
    feature_name varchar(50),
    metadata jsonb
);
$$

AFTER LOAD DO
$$
create index ix_load_rnc_sequence_features__accession on load_rnc_sequence_features(accession);
$$,
$$
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
)
$$,
$$
drop table load_rnc_sequence_features;
$$
;
