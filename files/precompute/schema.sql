DROP TABLE IF EXISTS precompute_urs_accession CASCADE;
CREATE TABLE IF NOT EXISTS precompute_urs_accession (
    id bigserial primary key,
    precompute_urs_id bigint not null,
    precompute_urs_taxid_id bigint not null,
    urs_taxid text not null,
    urs text not null,
    taxid int not null,
    is_active bool not null,
    last_release int not null,
    accession text not null,
    database text not null,
    description text not null,
    gene text,
    optional_id text,
    species text,
    common_name text,
    feature_name text,
    ncrna_class text,
    locus_tag text,
    organelle text,
    lineage text,
    so_rna_type text
);

DROP TABLE IF EXISTS precompute_urs_taxid CASCADE;
CREATE TABLE precompute_urs_taxid (
  id bigserial primary key,
  precompute_urs_id bigint not null,
  urs text not null,
  taxid int not null,
  urs_taxid text not null
);

DROP TABLE IF EXISTS precompute_urs CASCADE;
CREATE TABLE precompute_urs (
  id bigserial primary key,
  urs text not null
);

DROP TABLE IF EXISTS load_precomputed;
CREATE TABLE load_precomputed (
  id varchar(44) NOT NULL,
  upi varchar(26) NOT NULL,
  taxid int8 NULL,
  description varchar(500) NULL,
  short_description text NULL,
  rna_type varchar(500) NULL DEFAULT 'NULL'::character varying,
  has_coordinates bool NOT NULL DEFAULT false,
  databases text,
  is_active bool,
  last_release int4,
  so_rna_type text NOT NULL
);

DROP TABLE IF EXISTS load_qa_status;
CREATE TABLE load_qa_status (
  rna_id varchar(44) NOT NULL,
  upi varchar(26) NOT NULL,
  taxid int8 NOT NULL,
  has_issue bool not null,
  incomplete_sequence bool not null,
  possible_contamination bool not null,
  missing_rfam_match bool not null,
  from_repetitive_region bool not null,
  possible_orf bool not null,
  messages jsonb not null
);
