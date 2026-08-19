CREATE OR REPLACE FUNCTION upi.getupi(p_in_id bigint)
 RETURNS character
 LANGUAGE plpgsql
 IMMUTABLE SECURITY DEFINER
AS $function$
DECLARE
  v_upi RNACEN.xref.urs%TYPE := upper (to_hex ( p_in_id ) );
BEGIN

  ASSERT length (v_upi) <= 13, 'Can not generate upi: id ' || p_in_id || ' is too large';

  return substr ('URS0000000000', 1, 13 - length (v_upi) ) || v_upi;

END;
$function$
