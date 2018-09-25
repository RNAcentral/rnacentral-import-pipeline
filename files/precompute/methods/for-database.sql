drop table if exists upis_to_precompute;

CREATE TABLE upis_to_precompute as
select distinct
  upi
from xref
join rnc_database db on db.id = xref.dbid
where
  db.descr in (:dbs)
;

alter table upis_to_precompute
  add constraint fk_to_precompute__upi FOREIGN KEY (upi) REFERENCES rna(upi),
  add column id bigserial primary key
;
