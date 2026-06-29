CREATE OR REPLACE FUNCTION database.set_current_release(in_dbid bigint, in_release_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

      update rnc_database set current_release = in_release_id
         where id = in_dbid;
   end;

$function$

