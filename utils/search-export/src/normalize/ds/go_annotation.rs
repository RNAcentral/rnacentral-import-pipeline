use serde::{
    Deserialize,
    Serialize,
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GoAnnotation {
    pub id: usize,
    go_term_id: String,
    qualifier: String,
    go_name: String,
    assigned_by: String,
}
