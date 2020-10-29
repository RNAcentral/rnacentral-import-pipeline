DROP TABLE IF EXISTS load_locus;
CREATE TABLE load_locus (
  assembly_id text NOT NULL,
  locus_name text NOT NULL,
  locus_public_name text NOT NULL,
  chromosome text NOT NULL,
  strand text NOT NULL,
  locus_start int NOT NULL,
  locus_stop int NOT NULL,
  member_count int NOT NULL
  urs_taxid text NOT NULL,
  region_id int not null,
  member_type text NOT NULL
);

DROP TABLE IF EXISTS load_gene_status;
CREATE TABLE load_gene_status (
  assembly_id text NOT NULL,
  region_id int not null,
  urs_taxid text NOT NULL,
  cluster_status text NOT NULL
);

CREATE TABLE rnc_locus (
    id bigserial primary key,
    assembly_id text not null references ensembl_assembly(assembly_id),
    locus_name text not null,
    public_locus_name text NOT NULL UNIQUE,
    chromosome text not null,
    strand text not null,
    locus_start int not null,
    locus_stop int not null,
    member_count int not null,

    UNIQUE(assembly_id, locus_name)
);

CREATE TABLE rnc_locus_members (
    id bigserial primary key,
    urs_taxid text not null references rnc_rna_precomputed(id) ON DELETE CASCADE,
    region_id int not null UNIQUE references rnc_sequence_regions(id) ON DELETE CASCADE,
    locus_id bigint not null references rnc_locus(id) ON DELETE CASCADE,
    membership_status text NOT NULL
);
