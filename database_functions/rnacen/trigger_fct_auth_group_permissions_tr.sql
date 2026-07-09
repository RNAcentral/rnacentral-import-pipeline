CREATE OR REPLACE FUNCTION rnacen.trigger_fct_auth_group_permissions_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT auth_group_permissions_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

