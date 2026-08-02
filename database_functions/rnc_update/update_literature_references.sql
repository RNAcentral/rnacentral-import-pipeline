CREATE OR REPLACE FUNCTION rnc_update.update_literature_references()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    _context text;
    fcesig text;
    fcoid oid;
    sql_stmt text;
    -- Range-batch size (rows of the staging table per insert statement). Each
    -- statement queues one FK-check after-trigger event per inserted row
    -- (ck_rnc_reference_map__reference_id -> rnc_references.id, IMMEDIATE) and
    -- that queue is only drained at end of statement, so peak memory scales with
    -- batch size, NOT import size. A single unbatched insert OOMs in
    -- "AfterTriggerEvents".
    -- Smaller = safer/slower; this is the knob to turn down if it ever OOMs again.
    v_batch_size bigint := 10000000;
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

    -- Sort budget for the dedup/anti-join sorts below. SET LOCAL is txn-scoped so
    -- it doesn't leak to site connections. Kept modest; sorts spill to disk
    -- rather than getting a big budget.
    SET LOCAL work_mem = '256MB';

    -- update rnc_references table
    -- Only the distinct genuinely-new publications, so this stays small (bounded
    -- by the number of papers, not the number of xrefs) and needs no batching.
    sql_stmt := '
    INSERT INTO rnc_references(
      md5,
      location,
      authors,
      title,
      pmid,
      doi
    )
    SELECT
      in_md5,
      in_location,
      in_my_authors,
      in_title,
      in_pmid,
      in_doi
    FROM
      (
      WITH
        distinct_new_refs AS (
          SELECT DISTINCT
            md5,
            location,
            authors MY_AUTHORS,
            title,
            pmid,
            doi
          FROM load_rnc_references
         )
      SELECT
        l.md5 in_md5,
        l.location in_location,
        l.MY_AUTHORS in_my_authors,
        l.title in_title,
        l.pmid in_pmid,
        l.doi in_doi
      FROM rnc_references p RIGHT OUTER JOIN distinct_new_refs l ON (p.md5 = l.md5)
      WHERE p.md5 IS NULL
       ) alias4
';

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;

    -- update rnc_reference_map table
    --
    -- STEP 1 - stage the pairs to insert, once.
    -- A plain CREATE TABLE AS carries none of the FK/trigger/ON CONFLICT weight,
    -- so the single big join+sort here is safe. Two things happen in one pass:
    -- DISTINCT collapses the duplicate (accession, md5) rows the loader emits,
    -- and the anti-join drops pairs that already exist -- on a re-import that is
    -- almost all of them, and the unbatched insert used to pay a full ON CONFLICT
    -- index probe for every one. row_number() gives a dense integer key to slice
    -- on cheaply.
    RAISE NOTICE 'Building staging table load_ref_map';
    CREATE TEMP TABLE load_ref_map AS
    SELECT row_number() OVER () AS rn, d.accession, d.reference_id
    FROM (
      SELECT DISTINCT t3.accession, t4.id AS reference_id
      FROM rnacen.load_rnc_references t3
      JOIN rnacen.rnc_references t4 ON (t3.md5 = t4.md5)
      LEFT JOIN rnacen.rnc_reference_map m
        ON (m.accession = t3.accession AND m.reference_id = t4.id)
      WHERE m.accession IS NULL
    ) d;
    -- Tiny (one int column) index so each range slice is an indexed read, not a
    -- rescan of the whole staging table.
    CREATE INDEX ON load_ref_map(rn);
    ANALYZE load_ref_map;

    SELECT coalesce(max(rn), 0) INTO v_total FROM load_ref_map;
    RAISE NOTICE 'load_ref_map has % new pairs; inserting in batches of %',
                 v_total, v_batch_size;

    -- STEP 2 - range-batch the insert over disjoint rn slices.
    -- Each batch reads only its slice via the rn index; being a separate
    -- statement it fires + frees its FK-check queue when it completes, bounding
    -- peak memory to ~one batch. ON CONFLICT DO NOTHING stays as a safety net --
    -- the anti-join above is not atomic with the insert.
    lo := 1;
    WHILE lo <= v_total LOOP
        sql_stmt := format($q$
    insert into rnacen.rnc_reference_map as t1 (accession, reference_id)
    select accession, reference_id
      from load_ref_map
     where rn >= %s and rn < %s
        on conflict (accession, reference_id)
        do nothing
$q$, lo, lo + v_batch_size);

        -- EXPLAIN the first batch only, for the log (all batches share a plan).
        IF lo = 1 THEN
            explain_stmt := 'EXPLAIN (VERBOSE) ' || sql_stmt;
            EXECUTE explain_stmt into explain_result;
            RAISE NOTICE 'Batch plan: %', explain_result;
        END IF;

        RAISE NOTICE 'Inserting rn [%, %)', lo, lo + v_batch_size;
        EXECUTE sql_stmt;
        lo := lo + v_batch_size;
    END LOOP;

    -- TEMP table auto-drops at transaction end; explicit drop keeps peak disk low
    -- if anything runs after this in the same session.
    DROP TABLE load_ref_map;

    sql_stmt := '
truncate table load_rnc_references
';

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;


    RAISE NOTICE 'Literature references updated';

  END;

$function$
