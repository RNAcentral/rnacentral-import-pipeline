COPY (
SELECT
  json_build_object(
    'id', pre.urs_taxid,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rnc_rna_precomputed pre
JOIN rna
ON rna.urs = pre.urs
WHERE
  pre.rna_type = 'rRNA'
  and pre.is_active = true
  and pre.taxid is not null
  and (
    pre.databases ilike '%RDP%'
    or pre.databases ilike '%Ensembl%'
    or pre.databases ilike '%RefSeq%'
    or pre.databases ilike '%PDBe%'
    or pre.databases ilike '%FlyBase%'
    or pre.databases ilike '%MGI%'
    or pre.databases ilike '%PomBase%'
    or pre.databases ilike '%HGNC%'
    or pre.databases ilike '%SGD%'
    or pre.databases ilike '%RGD%'
    or pre.databases ilike '%TAIR%'
    or pre.databases ilike '%WormBase%'
  )
) TO STDOUT;
