CREATE OR REPLACE FUNCTION rnc_healthchecks.check_unique_upis()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE


    v_distinct_upi bigint;
    v_distinct_md5 bigint;
    v_all_rna_seq  bigint;

BEGIN

    SELECT count(DISTINCT urs) INTO v_distinct_upi FROM rna;
    SELECT count(DISTINCT md5) INTO v_distinct_md5 FROM rna;
    SELECT count(*)            INTO v_all_rna_seq  FROM rna;

    IF v_distinct_upi != v_distinct_md5 AND
       v_distinct_md5 != v_all_rna_seq THEN
      RAISE NOTICE 'not ok ... check_unique_upis';
    ELSE

      RAISE NOTICE 'ok ... check_unique_upis, found % RNA sequences', v_all_rna_seq;
    END IF;

  END;

$function$

