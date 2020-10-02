COPY (
  SELECT
    json_build_object(
      'urs', rna.upid,
      'length', rna.len,
      'md5', rna.md5
    )
  FROM rna
) TO STDOUT
