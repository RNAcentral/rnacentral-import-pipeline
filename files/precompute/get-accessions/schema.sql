CREATE TABLE IF NOT EXISTS urs_accession (
    id bigserial primary key,
    urs_taxid text not null,
    urs text not null,
    taxid int not null,
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
