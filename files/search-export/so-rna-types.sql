COPY (
select distinct
    coalesce(pre.assigned_so_rna_type, pre.so_rna_type)
FROM rnc_rna_precomputed pre
where
  coalesce(pre.assigned_so_rna_type, pre.so_rna_type) is not null
) TO STDOUT
