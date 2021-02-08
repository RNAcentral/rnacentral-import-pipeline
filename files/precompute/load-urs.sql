BEGIN TRANSACTION;

CREATE TEMP TABLE loaded_urs (
  urs_taxid TEXT NOT NULL
);

\copy loaded_urs from 'to-load.csv'

INSERT INTO precompute_urs (urs) (
SELECT DISTINCT
  split_part(urs_taxid, '_', 1)
from loaded_urs
);

/*
ALTER TABLE precompute_urs
  ADD CONSTRAINT un_precompute_urs__urs UNIQUE (urs),
  ADD CONSTRAINT fk_precompute_urs__urs FOREIGN KEY (urs) REFERENCES rna(upi)
;
*/

INSERT INTO precompute_urs_taxid (urs_taxid, precompute_urs_id, urs, taxid) (
SELECT 
  loaded.urs_taxid,
  urs.id,
  split_part(loaded.urs_taxid, '_', 1),
  split_part(loaded.urs_taxid, '_', 2)::int
FROM loaded_urs loaded
JOIN precompute_urs urs
ON
  urs.urs = split_part(loaded.urs_taxid, '_', 1)
);

ALTER TABLE precompute_urs_taxid
  ADD CONSTRAINT un_precompute_urs_taxid__urs_taxid UNIQUE (urs_taxid),
  ADD CONSTRAINT un_precompute_urs_taxid__urs__taxid UNIQUE (urs, taxid),
  ADD CONSTRAINT fk_precompute_urs_taxid__urs FOREIGN KEY (urs) REFERENCES rna(upi),
  ADD CONSTRAINT fk_precompute_urs_taxid__urs_id FOREIGN KEY(precompute_urs_id) REFERENCES precompute_urs(id)
;

COMMIT;
