use std::{
    error::Error,
    io,
    io::{
        BufRead,
        Write,
    },
    process::Output,
};

use tempfile::NamedTempFile;

/// This is a trait to extract the JSON encoded lines of something and produce a vector of all JSON
/// values parsed.
pub trait Jsonl {
    /// Parse and create a vector of each JSON encoded line. If any parsing step fails or there is
    /// an IO error then the whole process fails.
    fn jsonl(&self) -> Result<Vec<serde_json::Value>, Box<dyn Error>>;
}

impl Jsonl for Output {

    /// Parse the JSON encoded lines from stdout of some process.
    fn jsonl(&self) -> Result<Vec<serde_json::Value>, Box<dyn Error>> {
        let data = String::from_utf8_lossy(&self.stdout);
        let mut result = Vec::new();
        for line in data.lines() {
            result.push(serde_json::from_str(&line)?);
        }
        Ok(result)
    }
}

impl Jsonl for NamedTempFile {
    /// Parse the lines in the NamedTempFile to extract all JSON values. This does not seek to the
    /// start of the file, so that should be taken care of first.
    fn jsonl(&self) -> Result<Vec<serde_json::Value>, Box<dyn Error>> {
        let buf = io::BufReader::new(self);
        let mut data = Vec::new();
        for line in buf.lines() {
            let line = line?;
            data.push(serde_json::from_str(&line)?);
        }
        Ok(data)
    }
}

/// Build a NamedTempFile file with the lines of the given values.
pub fn temp_file_with(lines: Vec<&str>) -> io::Result<NamedTempFile> {
    let mut temp = NamedTempFile::new()?;
    for line in lines {
        writeln!(&mut temp, "{}", &line)?;
    }
    Ok(temp)
}
