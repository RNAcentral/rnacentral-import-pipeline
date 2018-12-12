DROP TABLE IF EXISTS :tablename;

CREATE TABLE :tablename as
select id, upi from rna;

alter table :tablename
  add constraint fk_to_precompute__upi FOREIGN KEY (upi) REFERENCES rna(upi),
  add primary key (id)
;
