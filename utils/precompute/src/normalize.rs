use serde::Deserialize;
use std::path::Path;

use anyhow::Result;

#[derive(Debug, Deserialize)]
pub struct Normalized {
    upi: String,
    taxid: Option<String>,
    length: usize,
    accessions: Vec<Accession>,
    coordinates: Vec<Coordinate>,
    deleted: Vec<bool>,
    previous: Previous,
    rfam_hits: Vec<RfamHit>,
    last_release: usize,
    r2dt_hits: Vec<R2dtHit>,
}

pub fn write(filename: &Path, output: &Path) -> Result<()> {
    let mut reader = rnc_utils::buf_reader(filename)?;
    let output = rnc_utils::buf_writer(output)?;
    let mut buf = String::new();
    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                let cleaned = buf.replace("\\\\", "\\");
                let raw: Raw = serde_json::from_str(&cleaned)?;
                let norm: Normalized = Normalized::try_from(raw)?;
                serde_json::to_writer(&mut output, &norm)?;
                writeln!(&mut output)?;
                buf.clear();
            }
        }
    }

    Ok(())
}
