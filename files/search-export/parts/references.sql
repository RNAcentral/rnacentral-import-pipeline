COPY (
  SELECT
    json_build_object(
      'id', xref.upi || '_' || xref.taxid,
      'authors', refs.authors,
      'journals', refs.location,
      'pub_titles', refs.title,
      'pub_ids', refs.id,
      'pubmed_ids', refs.pmid,
      'dois', refs.doi
    )
  FROM xref xref
  JOIN rnc_reference_map ref_map
  ON
    ref_map.accession = xref.ac
  JOIN rnc_references refs
  ON
    refs.id = ref_map.reference_id
) TO STDOUT

