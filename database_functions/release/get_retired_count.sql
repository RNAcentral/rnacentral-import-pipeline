CREATE OR REPLACE FUNCTION release.get_retired_count(in_dbid bigint, in_release bigint)
 RETURNS integer
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE


   v_count INTEGER := 0;

   v_previous RNACEN.rnc_release.id%TYPE;



BEGIN

   -- get id of previous database's release

    v_previous := release.get_previous_release(in_dbid, in_release);

   -- get the number of entries retired during the in_release

   select count(*) into v_count
      from RNACEN.xref
      where dbid = in_dbid and
            created < in_release and -- logically redundant but it optimizes the query

            last = v_previous;

   return v_count;

end;


/*
 * get_active_count
 */


$function$

