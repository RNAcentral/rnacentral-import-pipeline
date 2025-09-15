use bio::io::fasta;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Sequence<'a> {
    pub id: &'a str,
    pub description: Option<String>,
    pub sequence: &'a str,
}

impl<'a> From<Sequence<'a>> for fasta::Record {
    fn from(entry: Sequence<'a>) -> fasta::Record {
        fasta::Record::with_attrs(
            entry.id,
            entry.description.as_ref().map(|s| s.as_str()),
            entry.sequence.as_bytes(),
        )
    }
}
