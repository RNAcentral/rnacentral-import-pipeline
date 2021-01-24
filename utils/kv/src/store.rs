use std::{
    collections::HashMap,
    fs::File,
    io::{
        BufRead,
        BufReader,
    },
    path::{
        Path,
        PathBuf,
    },
    str,
    thread,
};

use serde_query::{
    DeserializeQuery,
    Query,
};

use serde_json::{
    Deserializer,
    Value,
};

use anyhow::{Result, Context};

use crossbeam_channel::{
    unbounded,
    Sender,
};

use rocksdb::{
    ColumnFamily,
    ColumnFamilyDescriptor,
    MergeOperands,
    Options,
    DB,
};

#[derive(DeserializeQuery)]
struct DocId {
    #[query(".id")]
    id: String,
}

pub struct Spec<'a> {
    path: &'a Path,
    allow_missing: bool,
    commit_size: usize,
    threads: usize,
}

impl<'a> Spec<'a> {
    pub fn new(path: &Path) -> Spec {
        Spec {
            path,
            allow_missing: false,
            commit_size: 1_000_000usize,
            threads: 4,
        }
    }

    pub fn set_allow_missing(&mut self, allow_missing: bool) -> () {
        self.allow_missing = allow_missing;
    }

    pub fn set_commit_size(&mut self, commit_size: usize) -> () {
        self.commit_size = commit_size;
    }
}

fn concat_merge(
    _new_key: &[u8],
    existing_val: Option<&[u8]>,
    operands: &mut MergeOperands,
) -> Option<Vec<u8>> {
    let mut result: Vec<u8> = Vec::with_capacity(operands.size_hint().0);
    existing_val.map(|v| {
        for e in v {
            result.push(*e)
        }
    });
    for op in operands {
        for e in op {
            result.push(*e)
        }
    }
    Some(result)
}

pub fn index(spec: &Spec, data_type: &str, filename: &Path) -> anyhow::Result<()> {
    let mut reader = rnc_utils::buf_reader(&filename)?;

    let mut db_opts = Options::default();
    db_opts.create_if_missing(true);
    db_opts.set_merge_operator("append operator", concat_merge, None);

    let mut families: Vec<String> = Vec::new();
    let mut store = match spec.path.exists() {
        true => {
            families.extend(DB::list_cf(&db_opts, spec.path)?);

            match families.len() {
                0 | 1 => {
                    DB::open(&db_opts, spec.path)
                },
                _ => {
                    let descriptors = families.clone().into_iter().map(|name| {
                        let mut cf_opts = Options::default();
                        cf_opts.set_merge_operator("append operator", concat_merge, None);
                        ColumnFamilyDescriptor::new(name, cf_opts)
                    });
                    DB::open_cf_descriptors(&db_opts, spec.path, descriptors)
                },
            }
        },
        false => DB::open(&db_opts, spec.path),
    }?;

    if !families.contains(&data_type.to_string()) {
        store.create_cf(data_type, &db_opts)?;
    }
    let family = store.cf_handle(data_type).unwrap();

    let mut buf = String::new();
    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                let line = buf.replace("\\\\", "\\");
                let id: DocId = serde_json::from_str::<Query<DocId>>(&line)?.into();
                store.merge_cf(&family, id.id.as_bytes(), line.as_bytes())?;
                buf.clear();
            },
        }
    }

    Ok(())
}

fn path_as_column_name(path: &Path) -> String {
    path.file_stem().unwrap().to_str().unwrap().to_string().replace("-", "_")
}

fn send_file_lines(path: &Path, sender: Sender<(String, DocId, String)>) {
    let file = File::open(path).unwrap();
    let data_type = path_as_column_name(&path);
    let mut reader = BufReader::new(file);
    let mut buf = String::new();
    loop {
        match reader.read_line(&mut buf).unwrap() {
            0 => break,
            _ => {
                let line = buf.replace("\\\\", "\\");
                let id: DocId = serde_json::from_str::<Query<DocId>>(&line).unwrap().into();
                sender.send((data_type.to_string(), id, line.to_string())).unwrap();
                buf.clear();
            },
        }
    }
}

pub fn index_files(spec: &Spec, filename: &Path) -> anyhow::Result<()> {
    let mut families: Vec<String> = Vec::new();
    let mut db_opts = Options::default();
    db_opts.create_if_missing(true);
    db_opts.set_merge_operator("append operator", concat_merge, None);

    let mut store = match spec.path.exists() {
        true => {
            families.extend(DB::list_cf(&db_opts, spec.path)?);
            families.sort();
            families.dedup();

            match families.len() {
                0 | 1 => DB::open(&db_opts, spec.path),
                _ => {
                    let descriptors = families.clone().into_iter().map(|name| {
                        let mut cf_opts = Options::default();
                        cf_opts.set_merge_operator("append operator", concat_merge, None);
                        ColumnFamilyDescriptor::new(name, cf_opts)
                    });
                    DB::open_cf_descriptors(&db_opts, spec.path, descriptors)
                },
            }
        },
        false => DB::open(&db_opts, spec.path),
    }?;

    let reader = rnc_utils::buf_reader(&filename)?;
    let (json_sender, json_reciever) = unbounded();
    for line in reader.lines() {
        let line = line?;
        let filename = line.trim_end();
        let path = PathBuf::from(filename);
        let column_name = path_as_column_name(&path);
        if !families.contains(&column_name) {
            store.create_cf(&column_name, &db_opts)?;
            families.push(column_name);
        }
        let sender = json_sender.clone();
        thread::spawn(move || send_file_lines(&path, sender));
    }

    let mut column_map: HashMap<String, &ColumnFamily> = HashMap::new();
    for name in families {
        column_map.insert(name.clone(), store.cf_handle(&name).unwrap());
    }

    for (data_type, id, data) in json_reciever {
        let column = column_map[&data_type];
        store.merge_cf(column, id.id.as_bytes(), data.as_bytes())?;
    }

    Ok(())
}

pub fn lookup(spec: &Spec, key_file: &Path, output: &Path) -> anyhow::Result<()> {
    let mut db_opts = Options::default();
    db_opts.set_merge_operator("append operator", concat_merge, None);
    let names = DB::list_cf(&db_opts, spec.path)
        .with_context(|| "Failed to get lookup names")?;
    let store = DB::open_cf_for_read_only(&db_opts, spec.path, &names, false)
        .with_context(|| "Failed to open read only db")?;

    let mut writer = rnc_utils::buf_writer(&output)
        .with_context(|| format!("Could not open {:?} for writing", &output))?;
    let mut keys = rnc_utils::buf_reader(&key_file)
        .with_context(|| format!("Could not open key file {:?}", &key_file))?;
    let mut buf = String::new();

    loop {
        match keys.read_line(&mut buf)? {
            0 => break,
            _ => {
                let mut data = HashMap::new();
                let trimmed = buf.trim_end();
                data.insert("id", serde_json::Value::String(trimmed.to_string()));
                let mut seen = false;
                let key = trimmed.as_bytes();
                for name in &names {
                    let default = serde_json::Value::Array(Vec::new());
                    let to_update = data.entry(&name).or_insert(default).as_array_mut().unwrap();
                    let cf = store.cf_handle(&name).unwrap();
                    match store.get_pinned_cf(cf, &key)? {
                        None => (),
                        Some(v) => {
                            seen = true;
                            let text = str::from_utf8(&v)?;
                            let values = Deserializer::from_str(text).into_iter::<Value>();
                            for value in values {
                                to_update.push(value?);
                            }
                        },
                    }
                }

                match (seen, spec.allow_missing) {
                    (true, _) => {
                        data.remove("__sled__default");
                    },
                    (false, true) => {
                        log::warn!("No data found for key {}", &buf);
                        continue;
                    },
                    (false, false) => {
                        return Err(anyhow::anyhow!("No data found for key {}", &buf));
                    },
                }

                serde_json::to_writer(&mut writer, &data)?;
                writeln!(&mut writer)?;
                buf.clear();
            },
        }
    }

    Ok(())
}
