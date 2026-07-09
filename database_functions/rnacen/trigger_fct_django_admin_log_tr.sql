CREATE OR REPLACE FUNCTION rnacen.trigger_fct_django_admin_log_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT django_admin_log_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

