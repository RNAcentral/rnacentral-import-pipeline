CREATE OR REPLACE FUNCTION rnacen.trigger_fct_django_comment_flags_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT django_comment_flags_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

