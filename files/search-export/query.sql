COPY (
  SELECT
    json_build_object(
        'rna_id', rna.upi || '_' || xref.taxid,
        'upi', rna.upi,
        'taxid', xref.taxid,
        'first_seen', array_agg(release1.timestamp),
        'last_seen', array_agg(release2.timestamp),
        'cross_references', array_agg(
            json_build_object(
                'name', acc."database",
                'external_id', acc.external_id,
                'optional_id', acc.optional_id,
                'accession', acc.accession,
                'non_coding_id', acc.non_coding_id,
                'parent_accession', acc.parent_ac || '.' || acc.seq_version
                )
            ),
        'description', array_agg(pre.description),
        'deleted', array_agg(xref.deleted),
        'length', array_agg(rna.len),
        'species', array_agg(acc.species),
        'organelles', array_agg(acc.organelle),
        'expert_dbs', array_agg(db.display_name),
        'rna_type', array_agg(pre.rna_type),
        'so_rna_type', array_agg(pre.so_rna_type),
        'product', array_agg(acc.product),
        'md5', array_agg(rna.md5),
        'authors', array_agg(refs.authors),
        'journals', array_agg(refs.location),
        'pub_titles', array_agg(refs.title),
        'pub_ids', array_agg(refs.id),
        'pubmed_ids', array_agg(refs.pmid),
        'dois', array_agg(refs.doi),
        'has_coordinates', array_agg(pre.has_coordinates),
        'rfam_family_names', array_agg(models.short_name),
        'rfam_ids', array_agg(hits.rfam_model_id),
        'rfam_clans', array_agg(models.rfam_clan_id),
        'qa_status', json_build_object(
            'has_issue', bool_or(qa.has_issue),
            'possible_contamination', bool_or(qa.possible_contamination),
            'incomplete_sequence', bool_or(qa.incomplete_sequence),
            'missing_rfam_match', bool_or(qa.missing_rfam_match)
        ),
        'tax_strings', array_agg(acc.classification),
        'functions', array_agg(acc.function),
        'genes', array_agg(acc.gene),
        'gene_synonyms', array_agg(acc.gene_synonym),
        'common_name', array_agg(acc.common_name),
        'notes', array_agg(acc.note),
        'locus_tags', array_agg(acc.locus_tag),
        'standard_names', array_agg(acc.standard_name),
        'products', array_agg(acc.product)
    )
FROM rna rna
join xref xref ON xref.upi = rna.upi
JOIN rnc_rna_precomputed pre ON xref.upi = pre.upi AND xref.taxid = pre.taxid
JOIN rnc_accessions acc ON xref.ac = acc.accession
JOIN rnc_database db ON xref.dbid = db.id
JOIN rnc_release release1 ON xref.created = release1.id
JOIN rnc_release release2 ON xref.last = release2.id
JOIN qa_status qa ON qa.rna_id = pre.id
LEFT JOIN rnc_reference_map ref_map ON ref_map.accession = acc.accession
LEFT JOIN rnc_references refs ON refs.id = ref_map.reference_id
LEFT JOIN rfam_model_hits hits ON xref.upi = hits.upi
LEFT JOIN rfam_models models ON hits.rfam_model_id = models.rfam_model_id
WHERE
    xref.deleted = 'N'
    AND rna.id BETWEEN :min AND :max
    AND pre.is_fragment = false
GROUP BY rna.upi, xref.taxid
) TO STDOUT
