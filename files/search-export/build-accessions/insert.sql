\timing

BEGIN TRANSACTION;

CREATE TEMP TABLE temp_current_urs_with_accession (
    id bigserial primary key,
    search_export_id int not null,
    urs_taxid text not null,
    urs text not null,
    taxid int not null,
    accession text not null
);

INSERT INTO temp_current_urs_with_accession
(search_export_id, urs_taxid, urs, taxid, accession)
(
SELECT
    todo.id,
    todo.urs_taxid,
    todo.urs,
    todo.taxid,
    xref.ac
FROM search_export_urs todo
JOIN :partition xref
ON
    xref.upi = todo.urs AND xref.taxid = todo.taxid
);

CREATE INDEX un_upis_accessions__accession ON temp_current_urs_with_accession(accession);

INSERT INTO search_export_accessions
(
    search_export_id,
    urs_taxid,
    urs,
    taxid,
    accession,
    lineage,
    common_name,
    database,
    external_id,
    function,
    gene,
    gene_synonym,
    locus_tag,
    non_coding_id,
    note,
    optional_id,
    organelle,
    parent_accession,
    product,
    species,
    standard_name
) (
SELECT
    todo.search_export_id,
    todo.urs_taxid,
    todo.urs,
    todo.taxid,
    todo.accession,
    acc.classification,
    acc.common_name,
    acc.database,
    acc.external_id,
    acc.function,
    acc.gene,
    acc.gene_synonym,
    acc.locus_tag,
    acc.non_coding_id,
    acc.note,
    acc.optional_id,
    acc.organelle,
    acc.parent_ac || '.' || acc.seq_version,
    acc.product,
    acc.species,
    acc.standard_name
FROM temp_current_urs_with_accession todo
JOIN rnc_accessions acc
ON
    acc.accession = todo.accession
WHERE
    acc.database = :'db_name'
);

COMMIT;
