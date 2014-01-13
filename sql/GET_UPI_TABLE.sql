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

create or replace function         get_upi_table(
   in_string in varchar2)
      return type_upi_table
   as
      v_string long := in_string || ',';

      v_pos number;

      v_upi_table type_upi_table := RNACEN.type_upi_table();

      begin
         loop
            v_pos := instr ( v_string, ',' );


            exit when ( nvl ( v_pos, 0 ) = 0 );

            v_upi_table.extend;

            v_upi_table ( v_upi_table.count ) := trim (
               substr ( v_string, 1, v_pos - 1 ) );

            v_string := substr( v_string, v_pos + 1);

         end loop;
         return (v_upi_table);
    end;
/
set define on
