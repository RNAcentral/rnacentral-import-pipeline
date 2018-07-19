COPY (
select
    json_build_object(
        'rna_id', pre.id,
        'rna_type',  pre.rna_type,
        'databases', pre.databases,
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
            'stop', mapping.stop
        )
    )
from xref
join rnc_rna_precomputed pre on pre.upi = xref.upi and pre.taxid = xref.taxid
LEFT JOIN rnc_coordinates coord ON xref.ac = coord.accession and coord.assembly_id = :assembly_id
LEFT JOIN rnc_genome_mapping mapping on mapping.rna_id = pre.id and mapping.assembly_id = :assembly_id
WHERE
    xref.taxid = :taxid
    and pre.is_active = true
    and xref.deleted = 'N'
    and (mapping.region_id is not null or coord.id is not null)
order by pre.id
) TO STDOUT
