COPY (
  SELECT
    json_build_object(
      'urs', rna.upi,
      'length', rna.len,
      'md5', rna.md5
    )
  FROM rna
) TO STDOUT
