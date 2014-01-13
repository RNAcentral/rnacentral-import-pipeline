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

set feedback off
set heading off
set termout off
set linesize 1000
set trimspool on
set verify off
spool &1..sql
prompt set define off
select decode( type||'-'||to_char(line,'fm99999'),
               'PACKAGE BODY-1', '/'||chr(10),
                null) ||
       decode(line,1,'create or replace ', '' ) ||
       text text
  from user_source
 where name = upper('&&1')
 order by type, line;
prompt /
prompt set define on
spool off
set feedback on
set heading on
set termout on
set linesize 100
