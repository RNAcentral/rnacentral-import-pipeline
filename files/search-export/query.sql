COPY (
SELECT
    json_build_object(
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
        'product', array_agg(acc.product),
        'md5', array_agg(rna.md5),
        'authors', array_agg(refs.authors),
        'journals', array_agg(refs.location),
        'pub_titles', array_agg(refs.title),
        'pub_ids', array_agg(refs.id),
        'pubmed_ids', array_agg(pubmed.ref_pubmed_id::varchar) || array_agg(refs.pmid),
        'dois', array_agg(pubmed.doi) || array_agg(refs.doi),
        'has_coordinates', array_agg(pre.has_coordinates),
        'rfam_family_names', array_agg(models.short_name),
        'rfam_ids', array_agg(hits.rfam_model_id),
        'rfam_clans', array_agg(models.rfam_clan_id),
        'rfam_status',
            case
                when cardinality((array_agg(pre.rfam_problems))) = 0 then '{}'
                else (array_agg(pre.rfam_problems))[1]::json
            end,
        'tax_strings', array_agg(acc.classification),
        'functions', array_agg(acc.function),
        'genes', array_agg(acc.gene),
        'gene_synonyms', array_agg(acc.gene_synonym),
        'common_name', array_agg(acc.common_name),
        'notes', array_agg(acc.note),
        'locus_tags', array_agg(acc.locus_tag),
        'standard_names', array_agg(acc.standard_name),
        'products', array_agg(acc.product),
        'go_annotations', array_agg(
            json_build_object(
                'go_term_id', anno.ontology_term_id,
                'qualifier', anno.qualifier,
                'go_name', ont.name,
                'assigned_by', anno.assigned_by
            )
        )
    )
FROM xref xref
JOIN rnc_accessions acc ON xref.ac = acc.accession
JOIN rnc_database db ON xref.dbid = db.id
JOIN rnc_release release1 ON xref.created = release1.id
JOIN rnc_release release2 ON xref.last = release2.id
JOIN rna rna ON xref.upi = rna.upi
JOIN rnc_rna_precomputed pre
ON
    xref.upi = pre.upi
    AND xref.taxid = pre.taxid
LEFT JOIN rnc_reference_map ref_map ON ref_map.accession = acc.accession
LEFT JOIN rnc_references refs ON refs.id = ref_map.reference_id
LEFT JOIN rfam_model_hits hits ON xref.upi = hits.upi
LEFT JOIN rfam_models models
ON
    hits.rfam_model_id = models.rfam_model_id
LEFT JOIN go_term_annotations anno ON anno.rna_id = pre.id
LEFT JOIN go_term_publication_map go_map
ON
    go_map.go_term_annotation_id = anno.go_term_annotation_id
LEFT JOIN ref_pubmed pubmed ON pubmed.ref_pubmed_id = go_map.ref_pubmed_id
LEFT JOIN ontology_terms ont
ON
    ont.ontology_term_id = anno.ontology_term_id
WHERE
  xref.deleted = 'N'
  AND rna.id BETWEEN :min AND :max
GROUP BY rna.upi, xref.taxid
) TO STDOUT
