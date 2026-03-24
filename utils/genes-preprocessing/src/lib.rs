use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
mod gene_preprocessing {
    use pyo3::prelude::*;
    use rayon::prelude::*;
    /// Formats the sum of two numbers as string.
    #[pyfunction]
    fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
        Ok((a + b).to_string())
    }

    #[pyfunction]
    fn count_overlap(
        exons_a: Vec<Vec<i64>>,
        exons_b: Vec<Vec<i64>>,
        threshold: f32,
    ) -> PyResult<i32> {
        let count = exons_a
            .iter()
            .zip(exons_b.iter())
            .filter_map(|(ea, eb)| {
                let start_a = ea[0];
                let end_a = ea[1];
                let start_b = eb[0];
                let end_b = eb[1];

                let length_a = end_a - start_a;
                let length_b = end_b - start_b;

                // Fix: overlap should be max(0, min_end - max_start)
                let overlap_start = start_a.max(start_b);
                let overlap_end = end_a.min(end_b);
                let overlap_length = (overlap_end - overlap_start).max(0);

                if overlap_length == 0 {
                    return None;
                }

                let overlap_ratio_a = overlap_length as f32 / length_a as f32;
                let overlap_ratio_b = overlap_length as f32 / length_b as f32;
                let min_ratio = overlap_ratio_a.min(overlap_ratio_b);

                if min_ratio >= threshold {
                    Some(())
                } else {
                    None
                }
            })
            .count() as i32;

        Ok(count)
    }

    #[pyfunction]
    fn count_overlap_batch(
        all_exons_a: Vec<Vec<Vec<i64>>>, // [pair_idx][exon_idx][start, end]
        all_exons_b: Vec<Vec<Vec<i64>>>,
        threshold: f32,
    ) -> PyResult<Vec<i32>> {
        let results = all_exons_a
            .par_iter()
            .zip(all_exons_b.par_iter())
            .map(|(exons_a, exons_b)| {
                // Your existing logic, just inlined
                exons_a
                    .iter()
                    .zip(exons_b.iter())
                    .filter(|(ea, eb)| {
                        let start_a = ea[0];
                        let end_a = ea[1];
                        let start_b = eb[0];
                        let end_b = eb[1];

                        let length_a = end_a - start_a;
                        let length_b = end_b - start_b;

                        let overlap_start = start_a.max(start_b);
                        let overlap_end = end_a.min(end_b);
                        let overlap_length = (overlap_end - overlap_start).max(0);

                        if overlap_length == 0 {
                            return false;
                        }

                        let overlap_ratio_a = overlap_length as f32 / length_a as f32;
                        let overlap_ratio_b = overlap_length as f32 / length_b as f32;
                        overlap_ratio_a.min(overlap_ratio_b) >= threshold
                    })
                    .count() as i32
            })
            .collect();

        Ok(results)
    }

    #[pyfunction]
    fn distance_to_agreement(exon_a_5p: Vec<i64>, exon_b_5p: Vec<i64>) -> PyResult<Vec<i64>> {
        Ok(vec![(exon_a_5p[0] - exon_b_5p[0]).abs(), (exon_a_5p[1] - exon_b_5p[1]).abs()])
        // Ok(exon_a_5p.iter().zip(exon_b_5p.iter()).map(|(xa, xb)|->
        // Vec<i64>{vec![(xa[0]-xb[0]).abs(), (xa[1]-xb[1]).abs()]}).collect())
    }

    #[pyfunction]
    fn distance_to_agreement_batch(
        exon_a_5p: Vec<Vec<i64>>,
        exon_b_5p: Vec<Vec<i64>>,
    ) -> PyResult<(Vec<i64>, Vec<i64>)> {
        let (diffs_0, diffs_1): (Vec<i64>, Vec<i64>) = exon_a_5p
            .par_iter()
            .zip(exon_b_5p.par_iter())
            .map(|(xa, xb)| ((xa[0] - xb[0]).abs(), (xa[1] - xb[1]).abs()))
            .unzip();
        Ok((diffs_0, diffs_1))
    }
}
