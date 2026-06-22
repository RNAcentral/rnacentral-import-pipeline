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

    -- update rnc_references table
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

    explain_stmt := 'EXPLAIN (VERBOSE) ' || sql_stmt;
    RAISE NOTICE 'Executing: %', explain_stmt;
    EXECUTE explain_stmt into explain_result;
    RAISE NOTICE '%', explain_result;

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;


    ----COMMIT;

    -- update rnc_reference_map table
    sql_stmt := '
    insert into rnc_reference_map as t1 (accession, reference_id)
    select t3.accession, t4.id 
      from load_rnc_references t3 
      join rnc_references t4 on (t3.md5 = t4.md5) 
        on conflict (accession, reference_id) 
        do nothing
';

    explain_stmt := 'EXPLAIN (VERBOSE) ' || sql_stmt;
    RAISE NOTICE 'Executing: %', explain_stmt;
    EXECUTE explain_stmt into explain_result;
    RAISE NOTICE '%', explain_result;

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;


    ----COMMIT;

    sql_stmt := '
truncate table load_rnc_references
';

    RAISE NOTICE 'Executing: %', sql_stmt;
    EXECUTE sql_stmt;


    RAISE NOTICE 'Literature references updated';

  END;

$function$

