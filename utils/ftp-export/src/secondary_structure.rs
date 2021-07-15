use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

#[derive(Debug, Serialize, Deserialize)]
pub struct SecondaryStructure {
    urs: String,
    sequence: String,
    secondary_structure: String,
}

pub fn write(raw: &Path, output: &Path) -> Result<()> {
    let mut reader = rnc_utils::buf_reader(raw)?;
    let mut writer = rnc_utils::buf_writer(output)?;
    let mut buf = String::new();
    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                let data: SecondaryStructure = serde_json::from_str(&buf)?;
                writeln!(&mut writer, ">{}", data.urs)?;
                writeln!(&mut writer, "{}", data.sequence)?;
                writeln!(&mut writer, "{}", data.secondary_structure)?;
                buf.clear();
            },
        }
    }

    writer.flush()?;
    Ok(())
}
