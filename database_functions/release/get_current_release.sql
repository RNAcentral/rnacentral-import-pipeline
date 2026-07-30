CREATE OR REPLACE FUNCTION release.get_current_release(in_dbid bigint)
 RETURNS bigint
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

      v_id RNACEN.rnc_release.id%TYPE;

BEGIN
      select current_release into STRICT  v_id from RNACEN.RNC_database
         where id = in_dbid;

      return v_id;
    exception
      when no_data_found then
         return v_id;
   end;


$function$

