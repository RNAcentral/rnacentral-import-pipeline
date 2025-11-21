use std::{
    io,
    path::{
        Path,
        PathBuf,
    },
    process::Output,
};

use serde_json::json;

use rnc_test_utils::{
    temp_file_with,
    Jsonl,
};

fn expand(id_file: &Path, json_file: &Path, output: &Path) -> io::Result<Output> {
    test_bin::get_test_bin("expand-urs").arg(id_file).arg(json_file).arg(output).output()
}

#[test]
fn simple_expanding_test() -> Result<(), Box<dyn std::error::Error>> {
    let id_file = temp_file_with(vec!["URS0000762A36_9606"])?;
    let json_file = temp_file_with(vec![r#"{"urs": "URS0000762A36", "value": "2"}"#])?;

    let output = PathBuf::from("-");
    let result = expand(id_file.path(), json_file.path(), &output)?;
    assert_eq!(String::from_utf8_lossy(&result.stderr), "");
    assert_eq!(
        result.jsonl()?,
        vec![json!({"id": "URS0000762A36_9606", "urs": "URS0000762A36", "value": "2"}),]
    );
    assert_eq!(result.status.success(), true);

    Ok(())
}
