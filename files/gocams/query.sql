-- URS_taxids to look up in GO-CAM. Directly-assayed annotations only
-- (ECO:0000314, "direct assay evidence"): GO-CAM models are built from curated
-- experimental evidence. Broaden the evidence_code filter if models turn up for
-- other evidence types.
SELECT DISTINCT urs_taxid
FROM go_term_annotations
WHERE evidence_code = 'ECO:0000314'
ORDER BY urs_taxid;
