CREATE OR REPLACE FUNCTION release.get_next_release(in_dbid bigint, in_release_id bigint)
 RETURNS bigint
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

      v_id RNACEN.rnc_release.id%TYPE;

BEGIN
      select min(id) into STRICT  v_id from RNACEN.rnc_release
         where dbid = in_dbid AND id > in_release_id;


      return v_id;
    exception
      when no_data_found then
         return v_id;
   end;


$function$

