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
    -- Number of hash buckets to split the upsert into. Peak per-statement memory
    -- is ~1/n_batches of the whole load, so bump this if a bigger import still OOMs.
    n_batches int := 16;
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

    -- Per-batch DISTINCT ON sort budget. SET LOCAL is txn-scoped so it doesn't
    -- leak to site connections. Kept modest (this VM is memory-tight); with
    -- batching each sort is only ~1/n_batches of the data anyway.
    SET LOCAL work_mem = '256MB';

    RAISE NOTICE 'Updating rnc_accessions in % batches', n_batches;

/*

    -- This index is not beneficial/used during the insert/update

    sql_stmt := '
create index if not exists load_rnc_accessions$accession on load_rnc_accessions(accession)
';
    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;
*/

    -- BATCHED UPSERT.
    -- A single INSERT ... ON CONFLICT over the whole load table (~228M rows)
    -- queues one FK-check after-trigger event per inserted row
    -- (rnc_accessions.rna_type -> ontology_terms) for the life of the statement,
    -- plus per-tuple ExecutorState allocations -> the process OOMs regardless of
    -- work_mem. Splitting into n_batches disjoint hash buckets means each
    -- statement fires + frees its FK checks and ExecutorState when it ends, so
    -- peak memory is ~1/n_batches and independent of how large the import is.
    --
    -- Bucketing on accession (the conflict arbiter key) keeps every row for a
    -- given accession in the SAME batch, so DISTINCT ON dedup and ON CONFLICT
    -- stay correct. abs(hashtextextended(...)) is a stable per-accession bucket.
    --
    -- ponytail: hash buckets mean each batch re-scans load_rnc_accessions (~n
    -- passes total). If that scan I/O becomes the bottleneck, switch to range
    -- batching on accession backed by a btree on load_rnc_accessions(accession)
    -- for a single ordered pass (and the DISTINCT ON sort disappears too).
    FOR i IN 0..(n_batches - 1) LOOP
        sql_stmt := format($q$
    insert into rnacen.rnc_accessions as t1 (
      accession, parent_ac, seq_version, feature_start, feature_end,
      feature_name, description, organelle, chromosome, function, gene,
      gene_synonym, inference, locus_tag, mol_type, ncRNA_class, note,
      product, standard_name, non_coding_id, database, external_id,
      optional_id, db_xref, rna_type, url
    )
    SELECT DISTINCT ON (t2.accession) accession, t2.parent_ac, t2.seq_version,
      t2.feature_start, t2.feature_end, t2.feature_name, t2.description,
      t2.organelle, t2.chromosome, t2.function, t2.gene, t2.gene_synonym,
      t2.inference, t2.locus_tag, t2.mol_type, t2.ncRNA_class, t2.note,
      t2.product, t2.standard_name, t2.non_coding_id,
      t2.database, t2.external_id, t2.optional_id, t2.db_xref, t2.so_term, t2.url
    FROM rnacen.load_rnc_accessions as t2
    WHERE mod(abs(hashtextextended(t2.accession, 0)), %s) = %s
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
    -- actually changed. Postgres does NOT support row-wise IS DISTINCT FROM
    -- between row constructors, so compare column-by-column and OR them.
    -- IS DISTINCT FROM is NULL-safe. This list MUST stay identical to the SET
    -- list above. Avoids rewriting every unchanged accession (dead
    -- tuples/WAL/bloat) on reloads where most rows are unchanged.
    WHERE t1.description    IS DISTINCT FROM excluded.description
       OR t1.organelle      IS DISTINCT FROM excluded.organelle
       OR t1.chromosome     IS DISTINCT FROM excluded.chromosome
       OR t1.function       IS DISTINCT FROM excluded.function
       OR t1.feature_name   IS DISTINCT FROM excluded.feature_name
       OR t1.gene           IS DISTINCT FROM excluded.gene
       OR t1.gene_synonym   IS DISTINCT FROM excluded.gene_synonym
       OR t1.inference      IS DISTINCT FROM excluded.inference
       OR t1.locus_tag      IS DISTINCT FROM excluded.locus_tag
       OR t1.mol_type       IS DISTINCT FROM excluded.mol_type
       OR t1.ncRNA_class    IS DISTINCT FROM excluded.ncRNA_class
       OR t1.note           IS DISTINCT FROM excluded.note
       OR t1.product        IS DISTINCT FROM excluded.product
       OR t1.standard_name  IS DISTINCT FROM excluded.standard_name
       OR t1.non_coding_id  IS DISTINCT FROM excluded.non_coding_id
       OR t1.database       IS DISTINCT FROM excluded.database
       OR t1.external_id    IS DISTINCT FROM excluded.external_id
       OR t1.optional_id    IS DISTINCT FROM excluded.optional_id
       OR t1.db_xref        IS DISTINCT FROM excluded.db_xref
       OR t1.rna_type       IS DISTINCT FROM excluded.rna_type
       OR t1.url            IS DISTINCT FROM excluded.url
$q$, n_batches, i);

        -- EXPLAIN the first batch only, for the log (all batches share a plan).
        IF i = 0 THEN
            explain_stmt := 'EXPLAIN (VERBOSE) ' || sql_stmt;
            EXECUTE explain_stmt into explain_result;
            RAISE NOTICE 'Batch plan: %', explain_result;
        END IF;

        RAISE NOTICE 'Executing batch % of %', i + 1, n_batches;
        EXECUTE sql_stmt;
    END LOOP;

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
