-- Copyright [2009-2014] EMBL-European Bioinformatics Institute
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

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
