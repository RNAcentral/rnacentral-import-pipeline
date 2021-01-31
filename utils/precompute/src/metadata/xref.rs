use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Deserialize, Serialize)]
pub struct Xref {
    #[serde(rename = "id")]
    pub urs_taxid: String,
    pub length: usize,
    pub deleted: String,
    pub last_release: usize,
    pub accession: String,
    pub ordering_index: usize,
}

impl Xref {
    pub fn is_active(&self) -> bool {
        self.deleted == "N"
    }
}
