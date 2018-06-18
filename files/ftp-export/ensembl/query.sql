COPY (
SELECT
    json_build_object(
        'rnacentral_id', pre.id,
        'description', max(pre.description),
        'sequence', COALESCE(rna.seq_short, rna.seq_long),
        'md5', max(rna.md5),
        'rna_type', max(pre.rna_type),
        'taxon_id', max(xref.taxid),
        'xrefs', array_agg(
            json_build_object(
                'database', db.display_name,
                'external_id', acc.external_id,
                'optional_id', acc.optional_id,
                'molecule_type', acc.mol_type
            )
        )
    )
FROM rna
JOIN xref ON xref.upi = rna.upi
JOIN rnc_rna_precomputed pre ON pre.upi = xref.upi AND pre.taxid = xref.taxid
JOIN rnc_database db ON db.id = xref.dbid
JOIN rnc_accessions acc
ON
    xref.ac = acc.accession
WHERE
    xref.deleted = 'N'
    AND rna.id BETWEEN {min_id} AND {max_id}
group by pre.id
) TO STDOUT
