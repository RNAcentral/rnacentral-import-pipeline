COPY (
select
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
from rna
left join pipeline_tracking_qa_scan analyzed
on
  analyzed.upi = rna.upi
  and analyzed.model_source = 'rfam'
where
  analyzed.upi is null
) TO STDOUT
