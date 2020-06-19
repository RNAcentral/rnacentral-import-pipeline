COPY (
select distinct
    pre.so_rna_type
FROM rnc_rna_precomputed pre
where
  pre.taxid is not null
  and pre.is_active = true
  and pre.so_rna_type is not null
) TO STDOUT
