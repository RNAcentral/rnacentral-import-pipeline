DROP TABLE IF EXISTS :tablename;
CREATE TABLE :tablename (
  id bigserial primary key,
  urs_taxid text not null,
  urs text not null,
  taxid int not null
);

create temp table loaded_urs (
  urs_taxid text not null
);

\copy loaded_urs from 'to-load.csv'

INSERT INTO :tablename (urs_taxid, urs, taxid) (
SELECT 
  urs_taxid,
  split_part(urs_taxid, '_', 1),
  split_part(urs_taxid, '_', 2)::int
from loaded_urs
);

ALTER TABLE :tablename 
  add constraint un_to_precompute__urs_taxid UNIQUE (urs_taxid),
  add constraint un_to_precompute__urs_taxid UNIQUE (urs, taxid),
  add constraint fk_to_precompute__urs FOREIGN KEY (urs) REFERENCES rna(upi)
;
