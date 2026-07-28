CREATE OR REPLACE FUNCTION rnc_healthchecks.check_expert_db_accessions()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_count bigint;
  dbs RECORD;

BEGIN

    FOR dbs IN (select descr from rnc_database where descr != 'ENA')
    LOOP
      SELECT count(*) INTO v_count FROM rnc_accessions WHERE database = dbs.descr AND external_id IS NULL;
      perform rnc_healthchecks.report_results(v_count, dbs.descr || ' expert accessions');
    END LOOP;

  END;

$function$

