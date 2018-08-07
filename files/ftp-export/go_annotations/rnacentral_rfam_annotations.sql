COPY (
select
    pre.id,
    pre.upi,
    pre.taxid,
    go_terms.ontology_term_id,
    string_agg('Rfam:' || hits.rfam_model_id, '|') as models
from rnc_rna_precomputed pre
join rfam_model_hits hits on hits.upi = pre.upi
join rfam_go_terms go_terms on hits.rfam_model_id = go_terms.rfam_model_id
where
    exists(select 1 from qa_status qa where qa.upi = pre.upi and qa.taxid = pre.taxid and qa.has_issue = false)
    and exists(select 1 from xref where xref.upi = pre.upi and xref.taxid = pre.taxid and xref.deleted = 'N')
group by pre.id, go_terms.ontology_term_id
) TO STDOUT
