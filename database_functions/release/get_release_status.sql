CREATE OR REPLACE FUNCTION release.get_release_status(in_release_id bigint)
 RETURNS character
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

      v_release_status RNACEN.rnc_release.status%TYPE;

BEGIN
      select status into STRICT  v_release_status
         from RNACEN.rnc_release where id = in_release_id;

      return v_release_status;
    exception
      when no_data_found then
         RAISE no_data_found USING MESSAGE = 'No release found with id: ' || in_release_id;
   end;



$function$

