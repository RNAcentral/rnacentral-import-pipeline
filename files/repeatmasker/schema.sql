-- One-off migration: create the RepeatMasker results table.
--
-- Run this once before the analyze/repeatmasker, precompute and search-export
-- pipelines so the repeatmasker search part and the precompute QA query can join
-- it even before RepeatMasker has populated any rows. The repeatmasker load ctl
-- also creates it (IF NOT EXISTS) at load time, so re-running is harmless.
CREATE TABLE IF NOT EXISTS repeatmasker_results (
  urs_taxid TEXT PRIMARY KEY,
  has_repeats bool,
  repeat_coverage float,
  repeat_count integer
);
