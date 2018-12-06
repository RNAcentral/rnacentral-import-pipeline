COPY (
select
  display_name
from rnc_database
where
  alive = 'Y'
  and num_sequences > 0
) TO STDOUT
