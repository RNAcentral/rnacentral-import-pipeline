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

set termout off
set heading off
set feedback off
set linesize 50
spool xtmpx.sql
select '@getcode ' || object_name
from user_objects
where object_type in ( 'PROCEDURE', 'FUNCTION', 'PACKAGE' )
/
spool off
spool install_all_code.sql
select '@' || object_name
from user_objects
where object_type in ( 'PROCEDURE', 'FUNCTION', 'PACKAGE' )
/
spool off
set heading on
set feedback on
set linesize 130
set termout on
@xtmpx.sql
