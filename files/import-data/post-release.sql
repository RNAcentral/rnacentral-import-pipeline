create table if not exists load_rnc_related_sequences (
  source_accession varchar(100) NOT NULL,
  target_accession varchar(100) NOT NULL,
  relationship_type text NOT NULL,
  methods text[]
);

create table if not exists load_rnc_sequence_features (
    accession varchar(100) NOT NULL,
    taxid int not null,
    start int not null,
    stop int not null,
    feature_name varchar(50),
    metadata jsonb
);

create table if not exists load_rnc_secondary_structure (
    rnc_accession_id varchar(100),
    secondary_structure text,
    md5 varchar(32)
);

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
);

INSERT INTO rnc_secondary_structure (
    rnc_accession_id,
    secondary_structure,
    md5
) (
select distinct
    rnc_accession_id,
    secondary_structure,
    md5
from load_rnc_secondary_structure
)
-- We can't include two ON CONFLICT statements so this will do extra inserts
ON CONFLICT (rnc_accession_id, md5) DO UPDATE SET
    rnc_accession_id = excluded.rnc_accession_id,
    secondary_structure = excluded.secondary_structure,
    md5 = excluded.md5
;

INSERT INTO rnc_related_sequences (
  source_urs_taxid,
  source_accession,
  target_urs_taxid,
  target_accession,
  relationship_type,
  methods
) (
select
  source.id,
  load.source_accession,
  target.id,
  load.target_accession,
  load.relationship_type::related_sequence_relationship,
  load.methods
from load_rnc_related_sequences load
join xref source_xref on source_xref.ac = load.source_accession
join rnc_rna_precomputed source
ON
  source.upi = source_xref.upi
  and source.taxid = source_xref.taxid
left join xref target_xref on target_xref.ac = load.target_accession
left join rnc_rna_precomputed target
on
  target.upi = target_xref.upi
  and target.taxid = target_xref.taxid
)
ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
SET
  methods = EXCLUDED.methods || rnc_related_sequences.methods
;

drop table load_rnc_related_sequences;
drop table load_rnc_secondary_structure;
drop table load_rnc_sequence_features;
