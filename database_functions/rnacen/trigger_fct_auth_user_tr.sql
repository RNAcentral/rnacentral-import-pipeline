CREATE OR REPLACE FUNCTION rnacen.trigger_fct_auth_user_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT auth_user_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

