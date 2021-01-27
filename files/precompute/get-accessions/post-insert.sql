CREATE INDEX ix_urs_accession__urs_taxid ON urs_accession(urs_taxid);

CLUSTER urs_accession USING ix_urs_accession__urs_taxid;
