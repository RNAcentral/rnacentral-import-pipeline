CREATE OR REPLACE FUNCTION rnc_test.run_tests()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    l_cntr bigint;

BEGIN

    -- truncate all tables in the beginning
    perform teardown;


    -- 17 test cases
    FOR l_cntr IN 1..17
    LOOP
      RAISE NOTICE 'Running test %', l_cntr;

      perform rnc_test.run_test(l_cntr);
    END LOOP;

  END;

$function$

