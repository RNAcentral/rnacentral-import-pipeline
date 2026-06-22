CREATE OR REPLACE FUNCTION rnacen.trigger_fct_auth_group_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT auth_group_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

