COPY (
  SELECT
    urs
  FROM pipeline_tracking_qa_scan
  where
    model_source = :'name'
) TO STDOUT WITH (FORMAT CSV)
