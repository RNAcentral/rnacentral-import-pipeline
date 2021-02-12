CREATE INDEX ix_search_accessions__urs_taxid ON search_export_accessions(urs_taxid);
CREATE INDEX ix_search_accessions__precompute_id ON search_export_accessions(search_export_id);

ALTER TABLE search_export_accessions
  ADD CONSTRAINT fk_search_accessions__search_id FOREIGN KEY (search_export_id) REFERENCES search_export_urs(id)
;
