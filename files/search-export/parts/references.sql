COPY (
  SELECT
    json_build_object(
      'id', xref.upi || '_' || xref.taxid,
      'authors', refs.authors,
      'journal', refs.location,
      'pub_title', refs.title,
      'pub_id', refs.id,
      'pubmed_id', refs.pmid,
      'doi', refs.doi
    )
  FROM xref xref
  JOIN rnc_reference_map ref_map
  ON
    ref_map.accession = xref.ac
  JOIN rnc_references refs
  ON
    refs.id = ref_map.reference_id
) TO STDOUT
