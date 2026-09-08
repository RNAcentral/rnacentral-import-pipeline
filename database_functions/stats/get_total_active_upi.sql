CREATE OR REPLACE FUNCTION stats.get_total_active_upi(in_release_id integer)
 RETURNS integer
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

   v_total_active INTEGER;

BEGIN

     select count(distinct(urs))
     into v_total_active
     from xref a
     where created <= in_release_id
     and last >= (SELECT max(id) from rnc_release b
                  where b.dbid = a.dbid and b.id <= in_release_id
                  and release_type = 'F');


     return v_total_active;

   END;

$function$
