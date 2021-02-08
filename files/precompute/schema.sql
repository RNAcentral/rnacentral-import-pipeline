DROP TABLE IF EXISTS precompute_urs;
CREATE TABLE precompute_urs (
  id bigserial primary key,
  urs text not null
);

DROP TABLE IF EXISTS precompute_urs_taxid;
CREATE TABLE precompute_urs_taxid (
  id bigserial primary key,
  precompute_urs_id int not null,
  urs text not null,
  taxid int not null,
  urs_taxid text not null
);

DROP TABLE IF EXISTS precompute_urs_accession;
CREATE TABLE IF NOT EXISTS precompute_urs_accession (
    id bigserial primary key,
    precompute_urs_id int not null,
    urs_taxid text not null,
    urs text not null,
    taxid int not null,
    is_active bool not null,
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
