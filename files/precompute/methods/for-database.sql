drop table if exists :tablename;

CREATE TABLE :tablename as
select distinct
  urs
from xref
join rnc_database db on db.id = xref.dbid
where
  db.descr in (:dbs)
;

alter table :tablename
  add constraint fk_to_precompute__upi FOREIGN KEY (urs) REFERENCES rna(urs),
  add column id bigserial primary key
;
