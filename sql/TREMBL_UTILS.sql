set define off

create or replace PACKAGE TREMBL_UTILS
AS
   FUNCTION get_sequence_chunk (p_upi_id VARCHAR2,slab INTEGER)
   RETURN   VARCHAR2
   ;
END;
/
create or replace PACKAGE BODY TREMBL_UTILS

AS
   FUNCTION get_sequence_chunk (p_upi_id VARCHAR2, slab INTEGER)
   RETURN VARCHAR2
   IS

      clobString VARCHAR2(4000);
      sub_start  NUMBER;
      chunksize  NUMBER:=4000;
   BEGIN
      sub_start := ((slab-1)*4000)+1;
      SELECT substr(seq_long,sub_start,chunksize)
      INTO   clobString
      FROM   rnacen.rna
      WHERE  upi = p_upi_id;
      RETURN clobString;
   END;
END TREMBL_UTILS;
/
set define on
