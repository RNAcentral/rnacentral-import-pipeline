ALTER TABLE :tablename
  ADD COLUMN id BIGSERIAL PRIMARY KEY,
  ADD CONSTRAINT fk_to_precompute__upi FOREIGN KEY (urs) REFERENCES rna(upi)
;

CREATE INDEX ix_to_precompute__urs ON :tablename(urs);
CREATE INDEX ix_to_precompute__taxid ON :tablename(taxid);
CREATE UNIQUE INDEX ix_to_precompute__urs_taxid ON :tablename(urs, taxid);
