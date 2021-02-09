use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Deserialize, Serialize)]
pub struct Basic {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    pub urs: String,
    pub taxid: usize,
    pub length: usize,
}
