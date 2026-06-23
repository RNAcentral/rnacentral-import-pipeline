CREATE OR REPLACE FUNCTION database.get_database_descr(in_dbid bigint)
 RETURNS character varying
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

      v_descr RNACEN.RNC_database.descr%TYPE;

BEGIN
      select descr into STRICT  v_descr from RNACEN.RNC_database
         where id = in_dbid;

      return v_descr;

    exception
      when no_data_found then
         return v_descr;
   end;


$function$

