CREATE OR REPLACE FUNCTION rnc_test.check_test12()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    l_count bigint;

    l_test_id varchar(30);

BEGIN

    l_test_id := 'check_test12_12,13,14,16';


    SELECT count(*) INTO l_count FROM rna;
    perform rnc_test.assertequals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref;
    perform rnc_test.assertequals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref

    WHERE upi = current_setting('RNC_TEST.FIRST_UPI')::varchar(13) AND deleted = 'Y' AND created = 1 AND version_I = 1;
    perform rnc_test.assertequals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi != current_setting('RNC_TEST.FIRST_UPI')::varchar(13) AND deleted = 'N' AND created = 2 AND version_I = 2;
    perform rnc_test.assertequals(l_test_id, 1, l_count);


  END;

$function$

