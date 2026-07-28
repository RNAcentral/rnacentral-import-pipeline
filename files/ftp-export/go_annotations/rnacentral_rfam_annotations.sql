COPY (
select
  json_build_object(
    'id', pre.urs_taxid,
    'upi', pre.urs,
    'taxid', pre.taxid,
    'ontology_term_id', go_terms.ontology_term_id,
    'models', string_agg(distinct 'Rfam:' || hits.rfam_model_id, '|')
  )
from rnc_rna_precomputed pre
join rfam_model_hits hits on hits.urs = pre.urs
join rfam_go_terms go_terms on hits.rfam_model_id = go_terms.rfam_model_id
where
    exists(select 1 from qa_status qa where qa.urs = pre.urs and qa.taxid = pre.taxid and qa.has_issue = false)
    -- and exists(select 1 from xref where xref.upi = pre.upi and xref.taxid = pre.taxid and xref.deleted = 'N')
    and pre.is_active = true
    and pre.databases != 'Rfam'
group by pre.urs_taxid, go_terms.ontology_term_id
) TO STDOUT
