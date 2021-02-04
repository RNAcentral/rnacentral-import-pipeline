COPY (
  SELECT
    json_build_object(
      'urs_taxid', todo.urs_taxid,
      'authors', refs.authors,
      'journals', refs.location,
      'pub_titles', refs.title,
      'pub_ids', refs.id,
      'pubmed_ids', refs.pmid,
      'dois', refs.doi
    )
  FROM search_export_accessions todo
  JOIN rnc_reference_map ref_map
  ON
    ref_map.accession = todo.accession
  JOIN rnc_references refs
  ON
    refs.id = ref_map.reference_id
  ORDER BY todo.search_export_id
) TO STDOUT

