CREATE OR REPLACE FUNCTION stats.get_active_xref(in_release_id integer)
 RETURNS integer
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

   v_active_xref INTEGER;

BEGIN

     select count(*)
     into v_active_xref
     from xref
     where created <= in_release_id
     and deleted = 'N';

   return v_active_xref;

  END;

$function$

