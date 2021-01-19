DROP TABLE IF EXISTS :tablename;
CREATE TABLE :tablename (
  id bigserial primary key,
  urs text NOT NULL
);

-- Select Skipped somehow
INSERT INTO :tablename (urs) (
SELECT DISTINCT
  xref.upi
FROM xref
LEFT JOIN rnc_rna_precomputed pre
ON
  pre.upi = xref.upi
  AND pre.taxid = xref.taxid
WHERE
  pre.id IS NULL
);

-- Select novel
INSERT INTO :tablename (urs) (
SELECT distinct
  xref.upi
FROM xref
LEFT JOIN rnc_rna_precomputed pre
ON
  pre.upi = xref.upi
  AND pre.taxid = xref.taxid
WHERE
  pre.last_release is null
);

-- Select outdataed
INSERT INTO :tablename (urs) (
SELECT distinct
  xref.upi
FROM xref
JOIN rnc_rna_precomputed pre
ON
  pre.upi = xref.upi
  and pre.taxid = xref.taxid
WHERE
  xref.last > pre.last_release
);

DELETE FROM basket
WHERE id IN
    (SELECT id
    FROM 
        (SELECT id,
         ROW_NUMBER() OVER( PARTITION BY fruit
        ORDER BY  id DESC ) AS row_num
        FROM basket ) t
        WHERE t.row_num > 1 );

ALTER TABLE :tablename
  ADD CONSTRAINT un_to_precompute_upi UNIQUE (upi),
  ADD CONSTRAINT fk_to_precompute__upi FOREIGN KEY (upi) REFERENCES rna(upi),
;
