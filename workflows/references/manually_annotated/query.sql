COPY(
  SELECT
    refs.pmid,
    xref.upi || '_' || xref.taxid
  FROM rnc_accessions acc
  JOIN xref ON xref.ac = acc.accession
  JOIN rnc_reference_map rmap ON rmap.accession = acc.accession
  JOIN rnc_references refs ON refs.id = rmap.reference_id
  WHERE
	  xref.dbid IN (14, 16, 18, 20, 23, 24, 27, 44, 48)
	  AND xref.deleted = 'N'
	  AND refs.pmid IS NOT NULL
) TO STDOUT (FORMAT CSV)
