CREATE OR REPLACE FUNCTION rnc_test.truncate_all_tables()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    perform teardown;
    EXECUTE 'TRUNCATE TABLE RNACEN.release_stats';
    EXECUTE 'TRUNCATE TABLE RNACEN.rnc_accessions';
    EXECUTE 'TRUNCATE TABLE RNACEN.rnc_references';
    RAISE NOTICE 'All tables truncated';
  end;

$function$

