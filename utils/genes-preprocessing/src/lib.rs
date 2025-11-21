use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
mod gene_preprocessing {
    use pyo3::prelude::*;
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
}
