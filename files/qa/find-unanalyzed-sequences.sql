COPY (
select
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
from rna
left join qa_analysis_status analyzed
on
  analyzed.urs = rna.upi
  and analyzed.qa_source = :'program'
  and analyzed.qa_source_version = :'version'
where
  analyzed.urs is null
) TO STDOUT
