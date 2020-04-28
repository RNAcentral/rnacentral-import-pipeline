DROP TABLE IF EXISTS :species_to_map;
CREATE TABLE :species_to_map AS (
  SELECT
    ensembl_url,
    assembly_id,
    taxid,
    division
  FROM ensembl_assembly
  WHERE
    division != 'EnsemblProtists'
    AND division != 'EnsemblFungi'
    AND EXISTS(SELECT 1
      FROM rnc_rna_precomputed pre
      WHERE
        pre.taxid = ensembl_assembly.taxid
        AND pre.has_coordinates IS false
        AND pre.is_active IS true
    )
);

CREATE INDEX ix__species_to_map__taxid ON :species_to_map(taxid);

DROP TABLE IF EXISTS :tablename;
CREATE TABLE :tablename (
  urs_taxid text PRIMARY KEY REFERENCES rnc_rna_precomputed(id),
  urs text NOT NULL REFERENCES rna(upi),
  taxid int NOT NULL references rnc_taxonomy(id),
  seq text NOT NULL,
  ensembl_url text NOT NULL,
  assembly_id text NOT NULL REFERENCES ensembl_assembly(assembly_id),
  division text,

  unique (urs_taxid, assembly_id)
);

CREATE INDEX ix__to_map__urs ON :tablename(urs);
CREATE INDEX ix__to_map__urs_taxid ON :tablename(urs_taxid);
CREATE INDEX ix__to_map__taxid ON :tablename(taxid);
CREATE INDEX ix__to_map__assembly_id ON :tablename(assembly_id);

INSERT INTO :tablename (
  urs,
  taxid,
  urs_taxid,
  seq,
  ensembl_url,
  assembly_id,
  division
) (
SELECT
  rna.upi,
  pre.taxid,
  pre.id,
  COALESCE(rna.seq_short, rna.seq_long),
  species.ensembl_url,
  species.assembly_id,
  species.division
FROM rna
JOIN rnc_rna_precomputed pre
ON
  pre.upi = rna.upi
JOIN :species_to_map species
ON 
  species.taxid = pre.taxid
WHERE
  rna.len BETWEEN :min_length AND :max_length
);

DELETE FROM :tablename to_map
USING pipeline_tracking_genome_mapping gm
WHERE
  gm.urs_taxid = to_map.urs_taxid
;

DELETE FROM :tablename to_map
USING rnc_sequence_regions coord
WHERE
  coord.urs_taxid = to_map.urs_taxid
;
