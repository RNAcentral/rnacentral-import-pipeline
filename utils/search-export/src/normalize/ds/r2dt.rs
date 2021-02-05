use serde::{
    Deserialize,
    Serialize,
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct R2dt {
    pub id: usize,
    secondary_structure_model: String,
    secondary_structure_source: String,
}

#[derive(Debug, PartialEq, Eq, Serialize)]
pub struct SecondaryStructure {
    model_name: String,
    model_source: String,
}
