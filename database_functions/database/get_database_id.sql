CREATE OR REPLACE FUNCTION database.get_database_id(in_descr character varying)
 RETURNS bigint
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

      v_id RNACEN.RNC_database.id%TYPE;

BEGIN
      select id into STRICT  v_id from RNACEN.RNC_database
         where UPPER(descr) = UPPER(in_descr);


      return v_id;
    exception
      when no_data_found then
         return v_id;
   end;


$function$

