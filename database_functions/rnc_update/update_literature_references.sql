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
BEGIN
    get diagnostics _context = pg_context;
    fcesig := substring(_context from 'function (.*?) line');
    fcoid := to_regprocedure(fcesig);
    raise notice 'executing function: % oid: %', fcesig, fcoid;
    execute 'set application_name = ''' || fcesig || '''';

    -- Sort budget for the DISTINCT/anti-join sort below. SET LOCAL is txn-scoped
    -- so it doesn't leak to site connections. Kept modest; the sort spills to
    -- disk rather than getting a big budget.
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
    -- The anti-join drops pairs that already exist -- on a re-import that is
    -- nearly all of them, and ON CONFLICT alone paid an index probe for every
    -- one. DISTINCT collapses the duplicate (accession, md5) rows the loader
    -- emits. ON CONFLICT DO NOTHING stays as a safety net: the anti-join is not
    -- atomic with the insert.
    sql_stmt := '
    insert into rnacen.rnc_reference_map as t1 (accession, reference_id)
    select distinct t3.accession, t4.id
      from rnacen.load_rnc_references t3
      join rnacen.rnc_references t4 on (t3.md5 = t4.md5)
      left join rnacen.rnc_reference_map m
        on (m.accession = t3.accession and m.reference_id = t4.id)
     where m.accession is null
        on conflict (accession, reference_id)
        do nothing
';

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;

    sql_stmt := '
truncate table load_rnc_references
';

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;


    RAISE NOTICE 'Literature references updated';

  END;

$function$
