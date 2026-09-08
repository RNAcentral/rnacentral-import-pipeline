CREATE OR REPLACE FUNCTION rnc_test.check_test8()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE


    l_count bigint;
    l_test_id varchar(20);

BEGIN

    l_test_id := 'check_test8';

    SELECT count(*) INTO l_count FROM rna;
    perform rnc_test.assertequals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    perform rnc_test.assertequals(l_test_id, 2, l_count);



    SELECT count(*) INTO l_count FROM xref
    WHERE urs = current_setting('RNC_TEST.FIRST_UPI')::varchar(13) AND deleted = 'Y' AND created = 1 AND LAST = 1 AND version = 1 AND taxid = 1 AND ac = 'id1';
    perform rnc_test.assertequals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE urs = current_setting('RNC_TEST.FIRST_UPI')::varchar(13) AND deleted = 'N' AND created = 2 AND LAST = 2 AND version = 2 AND taxid = 2 AND ac = 'id2';
    perform rnc_test.assertequals(l_test_id, 1, l_count);

  END;

$function$
