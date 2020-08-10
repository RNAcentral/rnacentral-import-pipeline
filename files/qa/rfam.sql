COPY (
select
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
from rna
left join pipeline_tracking_qa_scan analyzed
on
  analyzed.urs = rna.upi
  and analyzed.model_source = 'rfam'
where
  analyzed.urs is null
  and exists(select 1 from xref where xref.deleted = 'N' and rna.upi = xref.upi)
) TO STDOUT
