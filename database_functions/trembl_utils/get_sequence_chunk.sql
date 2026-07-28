CREATE OR REPLACE FUNCTION trembl_utils.get_sequence_chunk(p_upi_id text, slab numeric)
 RETURNS character varying
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE


      clobString varchar(4000);
      sub_start  bigint;
      chunksize  bigint:=4000;

BEGIN
      sub_start := ((slab-1)*4000)+1;
      SELECT substr(seq_long,sub_start,chunksize)
      INTO   clobString
      FROM   rnacen.rna
      WHERE  urs = p_upi_id;
      RETURN clobString;
   END;

$function$

