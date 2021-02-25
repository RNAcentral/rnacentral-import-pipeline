\timing

BEGIN;

CREATE INDEX ix_precompute_urs_accession__urs_taxid ON precompute_urs_accession(urs_taxid);
CREATE INDEX ix_precompute_urs_accession__precompute_id ON precompute_urs_accession(precompute_urs_id);

COMMIT;
