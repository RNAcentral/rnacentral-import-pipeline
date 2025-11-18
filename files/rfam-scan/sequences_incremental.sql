-- Select the next batch of sequences that will be scanned with Rfam
-- Selects from the pipeline_tracking_qa_scan table with a limit on the
-- number of sequences, sorted by the date

COPY (
  SELECT
    urs
  FROM pipeline_tracking_qa_scan
  where
    model_source = :'name'
  -- Try to handle semantic versioning: We allow the same major version (<=) and select things on an older minor version
  -- Have to use the regexp split to separate the two and then treat them as integers for comparison.
  -- This isn't fast, but seems to work ok for the R2DT table where I tested it.
  and (
    CAST((regexp_split_to_array(source_version, '\.'))[1] AS INT) <= CAST((regexp_split_to_array(:'min_version', '\.'))[1] AS INT)
    and
    CAST ((regexp_split_to_array(source_version, '\.'))[2] AS INT) < CAST((regexp_split_to_array(:'min_version', '\.'))[2] AS INT)
  )

  ORDER BY last_run
  LIMIT :'max_sequences'
) TO STDOUT WITH (FORMAT CSV)
