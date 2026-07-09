CREATE OR REPLACE FUNCTION release.set_release_status(in_release_id bigint, release_status character)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
      update RNACEN.rnc_release set status = release_status
         where id = in_release_id;
   end;



$function$

