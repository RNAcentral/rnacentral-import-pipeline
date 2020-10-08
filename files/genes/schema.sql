DROP TABLE IF EXISTS load_locus;
CREATE TABLE load_locus (
  taxid int not null,
  assembly_id text NOT NULL,
  locus_name text NOT NULL,
  chromosome text NOT NULL,
  strand text NOT NULL,
  locus_start int NOT NULL,
  locus_stop int NOT NULL,
  urs_taxid text NOT NULL,
  region_id int not null,
  is_representative bool NOT NULL
);
