create temp table temp_qa_to_scan AS
select
  rna.upi urs,
  COALESCE(rna.seq_short, rna.seq_long) as sequence
from rna
;

create index ix_temp_qa_to_scan on temp_qa_to_scan(urs);

delete from temp_qa_to_scan qa
using xref
where
  not exists(select 1 from xref where xref.deleted = 'N' and qa.urs = xref.upi)
;

delete from temp_qa_to_scan qa
using pipeline_tracking_qa_scan attempted
where
    attempted.urs = qa.urs
    and attempted.model_source = 'rfam'
;

COPY (
select
  json_build_object(
    'id', rna.urs,
    'sequence', rna.sequence,
  )
from temp_qa_to_scan rna
) TO STDOUT
