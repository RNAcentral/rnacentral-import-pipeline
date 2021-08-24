use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct GeneInfo {
    id: usize,
    name: String,
}

impl GeneInfo {
    pub fn id(&self) -> usize {
        self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }
}
