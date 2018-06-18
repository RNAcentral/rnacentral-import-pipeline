COPY (
select
    pre.id,
    pre.description,
    COALESCE(rna.seq_short, rna.seq_long)
from rnc_rna_precomputed pre
join rna on rna.upi = pre.upi
where
    pre.is_active = true
    and pre.taxid is null
) TO STDOUT
