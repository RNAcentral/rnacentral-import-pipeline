use serde::{
    Deserialize,
    Serialize,
};

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct Previous {
    pub id: usize,
    pub urs_taxid: String,
    upi: String,
    taxid: usize,
    databases: Option<String>,
    description: Option<String>,
    has_coordinates: Option<bool>,
    is_active: Option<bool>,
    last_release: Option<usize>,
    rna_type: Option<String>,
    short_description: Option<String>,
    so_rna_type: Option<String>,
}
