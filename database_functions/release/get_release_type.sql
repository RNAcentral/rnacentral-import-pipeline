CREATE OR REPLACE FUNCTION release.get_release_type(in_release_id bigint)
 RETURNS character
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE


      v_release_type RNACEN.rnc_release.release_type%TYPE;

BEGIN
      select release_type into STRICT  v_release_type
         from RNACEN.rnc_release where id = in_release_id;

      return v_release_type;
    exception
      when no_data_found then
         RAISE no_data_found USING MESSAGE = 'No release found with id: ' || in_release_id;
   end;


$function$

