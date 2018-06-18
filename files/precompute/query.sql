COPY (
select
    json_build_object(
        'upi', rna.upi,
        'taxid', xref.taxid,
        'accessions', array_agg(
            json_build_object(
                'description', acc.description,
                'gene', acc.gene,
                'optional_id', acc.optional_id,
                'pretty_database', db.display_name,
                'species', acc.species,
                'common_name', acc.common_name,
                'feature_name', acc.feature_name,
                'ncrna_class', acc.ncrna_class
            )
        ),
        'deleted', array_agg(xref.deleted),
        'xref_has_coordinates', array_agg(
            exists(
                select 1
                from rnc_coordinates coord
                where
                    coord.accession = xref.ac and xref.deleted = 'N'
            )),
        'rna_was_mapped', exists(
            select 1
            from rnc_genome_mapping mapping
            where
                mapping.upi = rna.upi and mapping.taxid = xref.taxid
        )
    )
FROM rna
JOIN xref
ON
    rna.upi = xref.upi
join rnc_database db
ON
    db.id = xref.dbid
LEFT JOIN rnc_accessions acc
ON
    xref.ac = acc.accession
    and xref.deleted = 'N'
where
  rna.id BETWEEN :min_id AND :max_id
GROUP BY rna.upi, xref.taxid
ORDER BY rna.upi
) TO STDOUT
