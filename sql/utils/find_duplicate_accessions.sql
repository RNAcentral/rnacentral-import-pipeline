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

WITH non_distinct_ac AS (
	SELECT
		count(*) as num,
		ac
	FROM load_rnacentral_all
	GROUP BY ac
	ORDER BY count(*) DESC
)
SELECT *
FROM load_rnacentral_all t, non_distinct_ac l
WHERE t.ac = l.ac AND l.num > 1
ORDER BY t.ac;
