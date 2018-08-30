COPY (
SELECT
    json_build_object(
        'rna_id', pre.id,
        'rna_type',  pre.rna_type,
        'databases', pre."databases",
        'known_coordinates', json_build_object(
            'region_id', coord.accession,
            'chromosome', coord.name,
            'strand', coord.strand,
            'start', coord.primary_start,
            'stop', coord.primary_end
        ),
        'mapped_coordinates', json_build_object(
            'region_id', mapping.region_id,
            'chromosome', mapping.chromosome,
            'strand', mapping.strand,
            'start', mapping."start",
            'stop', mapping.stop,
            'identity', mapping.identity
        )
    )
FROM rnc_rna_precomputed pre
JOIN xref
ON
    pre.upi = xref.upi AND pre.taxid = xref.taxid
JOIN ensembl_assembly assembly
ON
    assembly.taxid = xref.taxid
LEFT JOIN rnc_coordinates coord
ON
    xref.ac = coord.accession AND coord.assembly_id = :assembly_id
LEFT JOIN rnc_genome_mapping mapping
ON
    mapping.rna_id = pre.id AND mapping.assembly_id = :assembly_id
WHERE
    pre.is_active = true
    AND assembly.assembly_id = :assembly_id
    AND xref.deleted = 'N'
    AND (mapping.region_id IS NOT NULL OR coord.id IS NOT null)
-- ORDER BY pre.id
) TO STDOUT
