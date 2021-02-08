DROP TABLE IF EXISTS precompute_urs;
CREATE TABLE precompute_urs (
  id bigserial primary key,
  urs text not null,
);

DROP TABLE IF EXISTS precompute_urs_taxid;
CREATE TABLE precompute_urs_taxids (
  id bigserial primary key,
  precompute_urs_id int not null,
  urs text not null,
  taxid taxid not null,
  urs_taxid text not null
);
