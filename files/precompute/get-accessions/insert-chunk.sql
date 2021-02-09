BEGIN TRANSACTION;

CREATE TEMP TABLE temp_current_urs_with_accession (
    id bigserial primary key,
    precompute_urs_id int not null,
    precompute_urs_taxid_id int not null,
    urs_taxid text not null,
    urs text not null,
    taxid int not null,
    is_active bool not null,
    last_release int not null,
    accession text not null
);

INSERT INTO temp_current_urs_with_accession
(precompute_urs_id, precompute_urs_taxid_id, urs_taxid, urs, taxid, is_active, last_release, accession)
(
SELECT
    todo.precompute_urs_id,
    todo.id,
    todo.urs_taxid,
    todo.urs,
    todo.taxid,
    xref.deleted = 'N',
    xref.last_release,
    xref.ac
FROM precompute_urs_taxid todo
JOIN :partition xref
ON
    xref.upi = todo.urs AND xref.taxid = todo.taxid
);

CREATE INDEX un_upis_accessions__accession ON temp_current_urs_with_accession(accession);

INSERT INTO precompute_urs_accession
(
    precompute_urs_id,
    precompute_urs_taxid_id,
    urs_taxid,
    urs,
    taxid,
    is_active,
    last_release,
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
    todo.precompute_urs_id,
    todo.precompute_urs_taxid_id,
    todo.urs_taxid,
    todo.urs,
    todo.taxid,
    todo.is_active,
    todo.last_release,
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
