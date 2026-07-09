CREATE OR REPLACE FUNCTION rnc_update.update_rnc_accessions()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    _context text;
    fcesig text;
    fcoid oid;
    sql_stmt text;
-- the following could be a good idea
-- https://stackoverflow.com/questions/7682102/putting-explain-results-into-a-table
    explain_stmt text;
    explain_result text;
BEGIN
    get diagnostics _context = pg_context;
    fcesig := substring(_context from 'function (.*?) line');
    fcoid := to_regprocedure(fcesig);
    raise notice 'executing function: % oid: %', fcesig, fcoid;   
    execute 'set application_name = ''' || fcesig || '''';

    RAISE NOTICE 'Updating rnc_accessions';

/*

    -- This index is not beneficial/used during the insert/update

    sql_stmt := '
create index if not exists load_rnc_accessions$accession on load_rnc_accessions(accession)
';
    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;
*/

    sql_stmt := '
    insert into rnacen.rnc_accessions as t1 (
      accession, parent_ac, seq_version, feature_start, feature_end,
      feature_name, description,
      organelle, 
      chromosome, function, gene, gene_synonym, inference,
      locus_tag, mol_type, ncRNA_class, note,
      product, standard_name, non_coding_id,
      database, external_id, optional_id, db_xref, rna_type, url
    )
    SELECT DISTINCT ON (t2.accession) accession, t2.parent_ac, t2.seq_version, 
	  t2.feature_start, t2.feature_end, t2.feature_name, t2.description,
      t2.organelle, t2.chromosome, t2.function, t2.gene, t2.gene_synonym, 
	  t2.inference, t2.locus_tag, t2.mol_type, t2.ncRNA_class, t2.note,
      t2.product,  t2.standard_name, t2.non_coding_id,
      t2.database, t2.external_id, t2.optional_id, t2.db_xref, t2.so_term, t2.url
    FROM rnacen.load_rnc_accessions as t2
      on conflict (accession)
      do update
    SET
      description=excluded.description,
      organelle=excluded.organelle,
      chromosome=excluded.chromosome,
      function=excluded.function,
      feature_name=excluded.feature_name,
      gene=excluded.gene,
      gene_synonym=excluded.gene_synonym,
      inference=excluded.inference,
      locus_tag=excluded.locus_tag,
      mol_type=excluded.mol_type,
      ncRNA_class=excluded.ncRNA_class,
      note=excluded.note,
      product=excluded.product,
      standard_name=excluded.standard_name,
      non_coding_id=excluded.non_coding_id,
      database=excluded.database,
      external_id=excluded.external_id,
      optional_id=excluded.optional_id,
      db_xref=excluded.db_xref,
      rna_type=excluded.rna_type,
      url=excluded.url
    -- Skip no-op updates: only rewrite the row if at least one updated column
    -- actually changed. IS DISTINCT FROM is NULL-safe. This column list MUST stay
    -- identical to the SET list above. Avoids rewriting every unchanged accession
    -- (dead tuples/WAL/bloat) on reloads where most rows are unchanged.
    WHERE (
      t1.description, t1.organelle, t1.chromosome, t1.function, t1.feature_name,
      t1.gene, t1.gene_synonym, t1.inference, t1.locus_tag, t1.mol_type,
      t1.ncRNA_class, t1.note, t1.product, t1.standard_name, t1.non_coding_id,
      t1.database, t1.external_id, t1.optional_id, t1.db_xref, t1.rna_type, t1.url
    ) IS DISTINCT FROM (
      excluded.description, excluded.organelle, excluded.chromosome, excluded.function, excluded.feature_name,
      excluded.gene, excluded.gene_synonym, excluded.inference, excluded.locus_tag, excluded.mol_type,
      excluded.ncRNA_class, excluded.note, excluded.product, excluded.standard_name, excluded.non_coding_id,
      excluded.database, excluded.external_id, excluded.optional_id, excluded.db_xref, excluded.rna_type, excluded.url
    )
';


    explain_stmt := 'EXPLAIN (VERBOSE) ' || sql_stmt;
    RAISE NOTICE 'Executing: %', explain_stmt;
    EXECUTE explain_stmt into explain_result;
    RAISE NOTICE '%', explain_result;

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;

    --COMMIT;
    RAISE NOTICE 'rnc_accessions updated';

/*
    sql_stmt := '
drop index load_rnc_accessions$accession
';
    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;
*/

    sql_stmt := '
analyze rnacen.rnc_accessions
';
    RAISE NOTICE 'Executing: %', sql_stmt;
    execute sql_stmt;

    raise notice 'function: % oid: % completed', fcesig, fcoid;   
  END;

$function$

