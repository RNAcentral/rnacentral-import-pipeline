COPY (
select
    json_build_object(
        'upi', rna.upi,
        'taxid', xref.taxid,
        'length', rna.len,
        'accessions', array_agg(json_build_object(
                'is_active', xref.deleted = 'N',
                'description', acc.description,
                'gene', acc.gene,
                'optional_id', acc.optional_id,
                'pretty_database', db.display_name,
                'species', acc.species,
                'common_name', acc.common_name,
                'feature_name', acc.feature_name,
                'ncrna_class', acc.ncrna_class,
                'locus_tag', acc.locus_tag
            )
        ),
        'deleted', array_agg(distinct xref.deleted = 'Y'),
        'xref_has_coordinates', array_agg(exists(
            select 1
            from rnc_coordinates coord
            where
                coord.accession = xref.ac
        )),
        'rna_was_mapped', exists(
            select 1
            from rnc_genome_mapping mapping
            where
                mapping.upi = rna.upi
                and mapping.taxid = xref.taxid
        ),
        'previous', array_agg(row_to_json(prev.*))
    )
FROM rna
join xref on xref.upi = rna.upi
join rnc_accessions acc
on
	acc.accession = xref.ac
left join rnc_rna_precomputed prev
on
	prev.upi = rna.upi
	and prev.taxid = xref.taxid
join rnc_database db
ON
    db.id = xref.dbid
where
  rna.id BETWEEN :min AND :max
GROUP BY rna.upi, xref.taxid
ORDER BY rna.upi
) TO STDOUT
