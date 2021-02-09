use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct R2dtHit {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    model_id: usize,
    model_name: String,
    model_source: String,
    model_so_term: Option<String>,
    sequence_coverage: Option<f64>,
    model_coverage: Option<f64>,
    sequence_basepairs: Option<usize>,
    model_basepairs: Option<usize>,
}
