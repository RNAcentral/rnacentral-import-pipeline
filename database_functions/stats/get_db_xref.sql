CREATE OR REPLACE FUNCTION stats.get_db_xref(in_dbid integer)
 RETURNS integer
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE


  v_db_xref INTEGER;


BEGIN

    select count(*)
    into v_db_xref
    from xref
    where dbid = in_dbid
    and deleted = 'N';


  return v_db_xref;

  END;

$function$

