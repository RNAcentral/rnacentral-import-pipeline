CREATE OR REPLACE FUNCTION rnc_test.run_test(p_test_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    -- initialize rna, xref, and release tables
    perform setup;

    -- insert new data into the staging table

    -- +1 accounts for the first (default) entry
    perform rnc_test.import_staging_data( p_test_id + 1 );


    -- run update
    perform RNC_UPDATE.new_update;

    -- check results
    perform rnc_test.check_result(p_test_id);

    -- clean up
    perform teardown;


  END;

$function$

