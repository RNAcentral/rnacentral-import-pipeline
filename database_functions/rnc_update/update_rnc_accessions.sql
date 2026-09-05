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
    -- Range-batch size (rows of the deduped staging table per upsert statement).
    -- Each statement queues one FK-check after-trigger event per inserted row
    -- (rnc_accessions.rna_type -> ontology_terms) plus per-tuple ExecutorState for
    -- its whole life, so peak memory scales with batch size, NOT import size.
    -- Smaller = safer/slower; this is the knob to turn down if it ever OOMs again.
    v_batch_size bigint := 15000000;
    v_total bigint;
    lo bigint;
    explain_stmt text;
    explain_result text;
BEGIN
    get diagnostics _context = pg_context;
    fcesig := substring(_context from 'function (.*?) line');
    fcoid := to_regprocedure(fcesig);
    raise notice 'executing function: % oid: %', fcesig, fcoid;
    execute 'set application_name = ''' || fcesig || '''';

    -- Sort budget for the DISTINCT ON dedup below, raised from 256MB: at that size
    -- a ~110GB sort needs several merge passes, so the step is bound by temp I/O.
    -- populate_precompute.sql already runs this server at 2GB. Scoped to STEP 1 and
    -- given back before the upsert loop, whose trigger queue has OOM'd before.
    SET LOCAL work_mem = '2GB';
    -- CREATE INDEX reads maintenance_work_mem, not work_mem; unset it took the
    -- server default. 256MB is the value proven safe on rnc_rna_precomputed.
    SET LOCAL maintenance_work_mem = '256MB';

    -- STEP 1 - dedup once into a sorted staging table.
    -- A plain CREATE TABLE AS carries none of the FK/trigger/ON CONFLICT weight,
    -- so the single big sort here is safe (spills under work_mem, no per-row
    -- after-trigger queue). DISTINCT ON (accession) is the defensive dedup that
    -- must stay: it stops ON CONFLICT hitting the same accession twice in one
    -- command. Doing it ONCE here (instead of per hash bucket) means the whole
    -- 124 GB source is scanned + sorted exactly once. row_number() over the same
    -- accession ordering gives a dense integer key to slice on cheaply.
    RAISE NOTICE 'Building deduped staging table load_dedup';
    CREATE TEMP TABLE load_dedup AS
    SELECT row_number() OVER (ORDER BY accession) AS rn, d.*
    FROM (
      SELECT DISTINCT ON (accession)
        accession, parent_ac, seq_version, feature_start, feature_end,
        feature_name, description, organelle, chromosome, function, gene,
        gene_synonym, inference, locus_tag, mol_type, ncRNA_class, note,
        product, standard_name, non_coding_id, database, external_id,
        optional_id, db_xref, so_term, url
      FROM rnacen.load_rnc_accessions
      ORDER BY accession
    ) d;
    -- Tiny (one int column) index so each range slice is an indexed read, not a
    -- rescan of the ~110 GB staging table.
    CREATE INDEX ON load_dedup(rn);
    ANALYZE load_dedup;

    -- STEP 1 is done; give back the sort budget before the trigger-heavy loop.
    SET LOCAL work_mem = '256MB';

    SELECT max(rn) INTO v_total FROM load_dedup;
    RAISE NOTICE 'load_dedup has % deduped rows; upserting in batches of %',
                 v_total, v_batch_size;

    -- STEP 2 - range-batch the upsert over disjoint rn slices.
    -- Each batch reads only its slice via the rn index (no rescan, no re-sort);
    -- being a separate statement it fires + frees its FK-check queue and
    -- ExecutorState when it completes, bounding peak memory to ~one batch.
    -- rn slices are disjoint whole-accession sets (dedup already happened), so
    -- ON CONFLICT never sees the same accession twice within a batch.
    lo := 1;
    WHILE lo <= v_total LOOP
        sql_stmt := format($q$
    insert into rnacen.rnc_accessions as t1 (
      accession, parent_ac, seq_version, feature_start, feature_end,
      feature_name, description, organelle, chromosome, function, gene,
      gene_synonym, inference, locus_tag, mol_type, ncRNA_class, note,
      product, standard_name, non_coding_id, database, external_id,
      optional_id, db_xref, rna_type, url
    )
    SELECT accession, parent_ac, seq_version, feature_start, feature_end,
      feature_name, description, organelle, chromosome, function, gene,
      gene_synonym, inference, locus_tag, mol_type, ncRNA_class, note,
      product, standard_name, non_coding_id, database, external_id,
      optional_id, db_xref, so_term, url
    FROM load_dedup
    WHERE rn >= %s AND rn < %s
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
$q$, lo, lo + v_batch_size);

        -- EXPLAIN the first batch only, for the log (all batches share a plan).
        IF lo = 1 THEN
            explain_stmt := 'EXPLAIN (VERBOSE) ' || sql_stmt;
            EXECUTE explain_stmt into explain_result;
            RAISE NOTICE 'Batch plan: %', explain_result;
        END IF;

        RAISE NOTICE 'Upserting rn [%, %)', lo, lo + v_batch_size;
        EXECUTE sql_stmt;
        lo := lo + v_batch_size;
    END LOOP;

    RAISE NOTICE 'rnc_accessions updated';

    -- TEMP table auto-drops at transaction end; explicit drop keeps peak disk low
    -- if anything runs after this in the same session.
    DROP TABLE load_dedup;

    sql_stmt := '
analyze rnacen.rnc_accessions
';
    RAISE NOTICE 'Executing: %', sql_stmt;
    execute sql_stmt;

    raise notice 'function: % oid: % completed', fcesig, fcoid;
  END;

$function$
