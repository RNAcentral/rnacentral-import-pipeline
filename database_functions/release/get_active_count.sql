CREATE OR REPLACE FUNCTION release.get_active_count(in_dbid bigint, in_release bigint)
 RETURNS integer
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE


   v_count INTEGER := 0;


BEGIN

   -- get the number of active entries in the in_release

   select count(*) into v_count
      from RNACEN.xref
      where dbid = in_dbid 
        and created <= in_release 
        and ( deleted = 'N' or last >= in_release );

   return v_count;

end;

$function$

