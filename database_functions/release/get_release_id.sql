CREATE OR REPLACE FUNCTION release.get_release_id(in_dbid bigint, in_release_date timestamp without time zone)
 RETURNS bigint
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

      v_release_id RNACEN.rnc_release.id%TYPE;

BEGIN
      select id into STRICT  v_release_id
         from RNACEN.rnc_release where dbid = in_dbid and
            release_date = in_release_date;

      return v_release_id;
    exception

      when no_data_found then
         return v_release_id;
   end;


$function$

