COPY (
select
  json_build_object(
    'id', pre.id,
    'upi', pre.upi,
    'taxid', pre.taxid,
    'ontology_term_id', go_terms.ontology_term_id,
    'models', string_agg(distinct 'Rfam:' || hits.rfam_model_id, '|')
  )
from rnc_rna_precomputed pre
join rfam_model_hits hits on hits.upi = pre.upi
join rfam_go_terms go_terms on hits.rfam_model_id = go_terms.rfam_model_id
where
    exists(select 1 from qa_status qa where qa.upi = pre.upi and qa.taxid = pre.taxid and qa.has_issue = false)
    and exists(select 1 from xref where xref.upi = pre.upi and xref.taxid = pre.taxid and xref.deleted = 'N')
group by pre.id, go_terms.ontology_term_id
) TO STDOUT
