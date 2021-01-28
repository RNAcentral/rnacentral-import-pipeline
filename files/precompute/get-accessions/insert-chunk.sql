\timing

BEGIN TRANSACTION;

CREATE TEMP TABLE temp_current_urs_with_accession (
    id bigserial primary key,
    urs_taxid text not null,
    urs text not null,
    taxid int not null,
    is_active bool not null,
    accession text not null
);

INSERT INTO temp_current_urs_with_accession
(urs_taxid, urs, taxid, is_active, accession)
(
SELECT
    todo.urs_taxid,
    todo.urs,
    todo.taxid,
    xref.deleted = 'N',
    xref.ac
FROM upis_to_precompute todo
JOIN :partition xref
ON
    xref.upi = todo.urs AND xref.taxid = todo.taxid
);

CREATE INDEX un_upis_accessions__accession ON temp_current_urs_with_accession(accession);

INSERT INTO urs_accession
(
    urs_taxid,
    urs,
    taxid,
    accession,
    database,
    description,
    gene,
    optional_id,
    species,
    common_name,
    feature_name,
    ncrna_class,
    locus_tag,
    organelle,
    lineage,
    so_rna_type
) (
SELECT
    todo.urs_taxid,
    todo.urs,
    todo.taxid,
    todo.accession,
    acc.database,
    acc.description,
    acc.gene,
    acc.optional_id,
    acc.species,
    acc.common_name,
    acc.feature_name,
    acc.ncrna_class,
    acc.locus_tag,
    acc.organelle,
    acc.classification,
    acc.rna_type
FROM temp_current_urs_with_accession todo
JOIN rnc_accessions acc
ON
    acc.accession = todo.accession
WHERE
    acc.database = :'db_name'
);

COMMIT;