CREATE OR REPLACE FUNCTION rnc_test.check_test17()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE


    l_count bigint;
    l_test_id varchar(20);


BEGIN

    l_test_id := 'check_test17';

    SELECT count(*) INTO l_count FROM rna;
    perform rnc_test.assertequals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    perform rnc_test.assertequals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref

    WHERE urs = current_setting('RNC_TEST.FIRST_UPI')::varchar(13) AND deleted = 'Y' AND created = 1;

    perform rnc_test.assertequals(l_test_id, 1, l_count);

  END;

$function$

