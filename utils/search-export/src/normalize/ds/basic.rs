use serde::{
    Deserialize,
    Serialize,
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Basic {
    pub id: usize,
    pub length: usize,
    pub md5: String,
    pub urs_taxid: String,
}
