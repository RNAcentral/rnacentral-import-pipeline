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

select (select decode(extent_management,'LOCAL','*',' ') ||
               decode(segment_space_management,'AUTO','a ','m ')
	      from dba_tablespaces where tablespace_name = b.tablespace_name) || nvl(b.tablespace_name,
			 nvl(a.tablespace_name,'UNKOWN')) name,
	   kbytes_alloc kbytes,
	   kbytes_alloc-nvl(kbytes_free,0) used,
	   nvl(kbytes_free,0) free,
	   ((kbytes_alloc-nvl(kbytes_free,0))/
						  kbytes_alloc)*100 pct_used,
	   nvl(largest,0) largest,
	   nvl(kbytes_max,kbytes_alloc) Max_Size,
	   decode( kbytes_max, 0, 0, (kbytes_alloc/kbytes_max)*100) pct_max_used
from ( select sum(bytes)/1024 Kbytes_free,
			  max(bytes)/1024 largest,
			  tablespace_name
	   from  sys.dba_free_space
	   group by tablespace_name ) a,
     ( select sum(bytes)/1024 Kbytes_alloc,
			  sum(maxbytes)/1024 Kbytes_max,
			  tablespace_name
	   from sys.dba_data_files
	   group by tablespace_name
	   union all
      select sum(bytes)/1024 Kbytes_alloc,
			  sum(maxbytes)/1024 Kbytes_max,
			  tablespace_name
	   from sys.dba_temp_files
	   group by tablespace_name )b
where a.tablespace_name (+) = b.tablespace_name;
