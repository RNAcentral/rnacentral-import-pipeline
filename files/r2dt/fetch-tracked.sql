COPY(
  SELECT tracked.urs as urs
  FROM pipeline_tracking_traveler as tracked
) TO STDOUT CSV HEADER
