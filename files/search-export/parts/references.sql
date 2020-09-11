COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'authors', array_agg(refs.authors),
      'journals', array_agg(refs.location),
      'pub_titles', array_agg(refs.title),
      'pub_ids', array_agg(refs.id),
      'pubmed_ids', array_agg(refs.pmid),
      'dois', array_agg(refs.doi)
    )
  FROM :tablename todo
  JOIN xref xref
  ON
    xref.upi = todo.urs
  JOIN rnc_reference_map ref_map
  ON
    ref_map.accession = xref.ac
  JOIN rnc_references refs
  ON
    refs.id = ref_map.reference_id
  GROUP BY todo.id
) TO STDOUT

