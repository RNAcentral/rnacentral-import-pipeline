DROP TABLE IF EXISTS :tablename;
CREATE TABLE :tablename (
  id bigserial primary id,
  urs text not null
);

\copy :tablename from 'to-load.csv'

ALTER TABLE :tablename 
  add constraint un_to_precompute__urs UNIQUE (urs),
  add constraint fk_to_precompute__urs FOREIGN KEY (urs) REFERENCES rna(upi)
;
