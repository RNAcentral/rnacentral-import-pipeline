COPY (
select
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
from rna
left join rfam_analyzed_sequences analyzed
on
  analyzed.upi = rna.upi
where
  analyzed.upi is null
) TO STDOUT
