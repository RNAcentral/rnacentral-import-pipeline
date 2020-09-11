CREATE TABLE :tablename (
  id BIGSERIAL PRIMARY KEY,
  urs TEXT NOT NULL,
  taxid INT NOT NULL,
  urs_taxid TEXT NOT NULL
);

INSERT INTO :tablename (urs, urs_taxid, taxid) (
  SELECT
    pre.upi,
    pre.taxid,
    pre.id
  FROM xref
  WHERE
    pre.is_active = true
) ON CONFLICT DO NOTHING;

CREATE INDEX ix__search_todo_urs ON :tablename(urs);
CREATE INDEX ix__search_todo_taxid ON :tablename(taxid);
