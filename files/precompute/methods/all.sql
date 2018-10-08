drop table if exists upis_to_precompute;

CREATE TABLE upis_to_precompute as
select id, upi from rna;

alter table upis_to_precompute
  add constraint fk_to_precompute__upi FOREIGN KEY (upi) REFERENCES rna(upi),
  add primary key (id)
;
