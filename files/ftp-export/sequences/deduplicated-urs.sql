COPY (
  SELECT urs
  FROM rnc_sequence_cluster
  WHERE in_deduplicated_set
) TO STDOUT
