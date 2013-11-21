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
