use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct RfamHit {
    pub id: usize,
    pub urs_id: usize,
    urs_taxid: String,
    rfam_hit_id: usize,
    model: String,
    model_rna_type: Option<String>,
    model_domain: Option<String>,
    model_name: String,
    model_long_name: String,
    model_completeness: f64,
    model_start: usize,
    model_stop: usize,
    sequence_completeness: f64,
    sequence_start: usize,
    sequence_stop: usize,
}
