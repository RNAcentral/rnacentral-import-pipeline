CREATE OR REPLACE FUNCTION rnc_update.verify_xref_id_not_null()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_count bigint;

BEGIN
  -- update all null id
  UPDATE xref SET id = nextval('xref_pk_seq') WHERE id IS NULL;

  --COMMIT;

END;

$function$

