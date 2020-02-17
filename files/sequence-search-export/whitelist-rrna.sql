COPY (
SELECT
  json_build_object(
    'id', pre.id,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rnc_rna_precomputed pre
JOIN rna
ON rna.upi = pre.upi
WHERE
  pre.rna_type = 'rRNA'
  and pre.is_active = true
  and pre.taxid is not null
  and (
    pre.databases like '%RDP%'
    or pre.databases like '%Ensembl%'
    or pre.databases like '%RefSeq%'
    or pre.databases like '%PDBe%'
    or pre.databases like '%FlyBase%'
    or pre.databases like '%MGI%'
    or pre.databases like '%PomBase%'
    or pre.databases like '%HGNC%'
    or pre.databases like '%SGD%'
    or pre.databases like '%RGD%'
    or pre.databases like '%TAIR%'
    or pre.databases like '%WormBase%'
  )
) TO STDOUT;
