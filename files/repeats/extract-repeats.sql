SELECT sr.name, rf.seq_region_start - 1, rf.seq_region_end, logic_name
FROM seq_region sr
INNER JOIN repeat_feature rf USING (seq_region_id)
INNER JOIN analysis USING (analysis_id)
ORDER BY sr.name, rf.seq_region_start;
