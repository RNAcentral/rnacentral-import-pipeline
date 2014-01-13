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

create or replace procedure add_column is

   CURSOR c_names IS
   SELECT table_name FROM user_tables
   WHERE table_name like 'LOAD_%';

begin
   FOR v_names IN c_names LOOP

       if v_names.table_name != 'LOAD_ENSEMBL_ARMADILLO' THEN

       execute immediate


       'alter table '||v_names.table_name ||' add TAXID VARCHAR2(30)';

       end if;

   END LOOP;

end add_column;
/
set define on
