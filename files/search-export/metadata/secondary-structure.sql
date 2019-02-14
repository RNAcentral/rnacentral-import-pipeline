COPY (
SELECT
  json_build_object(
    'rna_id', pre.id,
    'secondary', json_build_object(
      'has_secondary', true
    )
  )
FROM rnc_secondary_structure_layout as layout
join rnc_rna_precomputed pre on pre.upi = layout.urs
WHERE
  pre.taxid is not null
UNION
select
  json_build_object(
    'rna_id', xref.upi || '_' || xref.taxid
    'secondary', json_build_object(
      'has_secondary', true
    )
  )
from rnc_secondary_structure as secondary
join xref on xref.ac = secondary.accession
where
  xref.deleted = 'N'
) TO STDOUT
