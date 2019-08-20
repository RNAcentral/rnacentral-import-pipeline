-- Update layout to hide any sequences with too small sequence coverage
UPDATE rnc_secondary_structure_layout layout
SET
  should_show = false
FROM rnc_secondary_structure_layout_models models
WHERE
  models.id = layout.model_id
  AND models.model_name ilike 'RF%'
  AND layout.sequence_coverage < 0.5
;

-- Hide all bad SSU
UPDATE rnc_secondary_structure_layout
SET
  should_show = false
FROM rnc_secondary_structure_layout_models models
WHERE
  models.id = layout.model_id
  AND models.so_term_id = 'SO:0000650'
  AND layout.overlap_count >= 14
;

-- Hide all introns
UPDATE rnc_secondary_structure_layout
SET
  should_show = false
FROM rnc_secondary_structure_layout_models models
WHERE
  models.id = layout.model_id
  AND models.so_term_id = 'SO:0000587'
;

-- Hide bad 5S
UPDATE rnc_secondary_structure_layout
SET
  should_show = false
FROM rnc_secondary_structure_layout_models models
WHERE
  models.id = layout.model_id
  AND models.so_term_id = 'SO:0000652'
  AND layout.overlap_count >= 3
;

-- Hide all Rfam families with many overlaps
UPDATE rnc_secondary_structure_layout
SET
  should_show = false
FROM rnc_secondary_structure_layout_models models
WHERE
  models.id = layout.model_id
  AND model.model_name ilike 'RF%'
  AND layout.overlap_count >= 15
;
