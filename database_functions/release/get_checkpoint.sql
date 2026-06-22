CREATE OR REPLACE FUNCTION release.get_checkpoint()
 RETURNS bigint
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

     v_checkpoint RNACEN.rnc_release.id%TYPE;

BEGIN
     select max(id) into v_checkpoint
        from RNACEN.rnc_release where status = 'D';

     return v_checkpoint;
  end;



/*
 * get_retired_count
 */


$function$

