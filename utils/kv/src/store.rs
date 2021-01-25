use std::{
    collections::HashMap,
    error::Error,
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

use anyhow::{
    Context,
    Result,
};

// use crossbeam_channel::{
//     unbounded,
//     Sender,
// };

use rkv::{
    backend::{
        BackendEnvironmentBuilder,
        Lmdb,
        LmdbDatabase,
        LmdbEnvironment,
    },
    Manager,
    MultiStore,
    Rkv,
    SingleStore,
    StoreOptions,
    Value,
};

#[derive(Debug, DeserializeQuery)]
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

pub fn index(spec: &Spec, data_type: &str, filename: &Path) -> Result<(), Box<dyn Error>> {
    let mut reader = rnc_utils::buf_reader(&filename)?;

    let mut manager = Manager::<LmdbEnvironment>::singleton().write()?;
    let mut builder = Rkv::environment_builder::<Lmdb>();
    builder.set_max_dbs(20);
    builder.set_make_dir_if_needed(true);
    builder.set_map_size(200000000000);
    let arc = manager
        .get_or_create_from_builder(spec.path, builder, Rkv::from_builder::<Lmdb>)
        .with_context(|| "Failed to create arc")?;
    let env = arc.read().unwrap();
    let store = env
        .open_single(data_type, StoreOptions::create())
        .with_context(|| format!("Failed to open store {:?}", data_type))?;

    let mut writer = env.write().with_context(|| "Failed to create writer")?;
    let mut buf = String::new();
    let mut count = 0usize;
    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                let line = buf.replace("\\\\", "\\");
                let line = line.trim_end();
                let id: DocId = serde_json::from_str::<Query<DocId>>(&line)?.into();
                let value = serde_json::from_str(&line)?;
                let existing = store.get(&writer, &id.id)?;
                let to_add = match existing {
                    None => {
                        let value = vec![value];
                        serde_json::to_string(&value)?
                    },
                    Some(v) => match v {
                        Value::Json(e) => {
                            let mut current: serde_json::Value = serde_json::from_str(&e)?;
                            let current = current.as_array_mut().unwrap();
                            current.push(value);
                            serde_json::to_string(&current)?
                        },
                        _ => panic!("Invalid existing data, should never happen"),
                    },
                };

                store.put(&mut writer, &id.id, &Value::Json(&to_add))?;
                buf.clear();
                count += 1;
                if count % spec.commit_size == 0 {
                    println!("Commiting count: {}", count);
                    writer
                        .commit()
                        .with_context(|| format!("Failed to write data from {:?}", filename))?;
                    writer = env.write().with_context(|| "Could not create for new chunk")?;
                }
            },
        }
    }
    writer.commit().with_context(|| format!("Failed to write data from {:?}", filename))?;

    Ok(())
}

// fn path_as_column_name(path: &Path) -> String {
//     path.file_stem().unwrap().to_str().unwrap().to_string().replace("-", "_")
// }

// fn send_file_lines(path: &Path, sender: Sender<(String, DocId, String)>) {
//     let file = File::open(path).unwrap();
//     let data_type = path_as_column_name(&path);
//     let mut reader = BufReader::new(file);
//     let mut buf = String::new();
//     loop {
//         match reader.read_line(&mut buf).unwrap() {
//             0 => break,
//             _ => {
//                 let line = buf.replace("\\\\", "\\");
//                 let id: DocId =
// serde_json::from_str::<Query<DocId>>(&line).unwrap().into();
// sender.send((data_type.to_string(), id, line.to_string())).unwrap();
// buf.clear();             },
//         }
//     }
// }

// pub fn index_files(spec: &Spec, filename: &Path) -> Result<(), Box<dyn Error>> {
//     let mut manager = Manager::<LmdbEnvironment>::singleton().write()?;
//     let arc = manager.get_or_create(spec.path, Rkv::new::<Lmdb>)?;
//     let env = arc.read();

//     // let stores = HashMap::<String, Multi>
//     let reader = rnc_utils::buf_reader(&filename)?;
//     let (json_sender, json_reciever) = unbounded();
//     for line in reader.lines() {
//         let line = line?;
//         let filename = line.trim_end();
//         let path = PathBuf::from(filename);
//         let sender = json_sender.clone();
//         let store = env
//             .open_multi(&data_type, StoreOptions::create())
//             .with_context(|| format!("Failed to open store {:?}", data_type))?;
//         stores.insert(name, store);
//         thread::spawn(move || send_file_lines(&path, sender));
//     }

//     let writer = env.write().with_context(|| "Failed to create writer")?;

//     for (data_type, id, data) in json_reciever {
//         let store = stores.get(&data_type)?;
//         store
//             .put(&mut writer, id.id, &Value::Json(data))
//             .with_context(|| format!("Failed to write data for {:}/{:?}", &data_type,
// &id))?;     }

//     Ok(())
// }

pub fn lookup(spec: &Spec, key_file: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
    let mut writer = rnc_utils::buf_writer(&output)
        .with_context(|| format!("Could not open {:?} for writing", &output))?;
    let mut keys = rnc_utils::buf_reader(&key_file)
        .with_context(|| format!("Could not open key file {:?}", &key_file))?;

    // let mut manager = Manager::<LmdbEnvironment>::singleton().write()?;
    // let mut builder = Rkv::environment_builder::<Lmdb>();
    // builder.set_max_dbs(20);
    // builder.set_make_dir_if_needed(true);
    // builder.set_map_size(200000000000);
    // println!("Created BUILDER");
    // let arc = manager
    //     .get_or_create_from_builder(spec.path, builder, Rkv::from_builder::<Lmdb>)
    //     .with_context(|| "Failed to create arc")?;
    // let env = arc.read().unwrap();

    let mut manager = Manager::<LmdbEnvironment>::singleton().write().unwrap();
    let created_arc = manager.get_or_create(spec.path, Rkv::new::<Lmdb>).unwrap();
    let env = created_arc.read().unwrap();

    let names = env.get_dbs()?.iter().filter_map(|db| db.to_owned()).collect::<Vec<String>>();

    let mut stores: Vec<(String, SingleStore<LmdbDatabase>)> = Vec::new();
    for name in names {
        let store = env
            .open_single(name.as_str(), StoreOptions::default())
            .with_context(|| format!("Failed to open store {:?}", &name))?;
        stores.push((name.to_owned(), store));
    }

    let reader = env.read().with_context(|| "Failed to create writer")?;

    let mut buf = String::new();
    loop {
        match keys.read_line(&mut buf)? {
            0 => break,
            _ => {
                let mut data = HashMap::new();
                let id = buf.trim_end().to_string();
                data.insert("id", serde_json::Value::String(id.to_string()));
                let mut seen = false;
                for (name, store) in &stores {
                    let default = serde_json::Value::Array(Vec::new());
                    let to_update = data.entry(&name).or_insert(default).as_array_mut().unwrap();
                    let current = store
                        .get(&reader, &id)
                        .with_context(|| format!("Failed to get results for {:?}", &id))?;

                    match current {
                        None => (),
                        Some(v) => match v {
                            Value::Json(raw) => {
                                seen = true;
                                let parsed = serde_json::from_str(&raw)?;
                                to_update.push(parsed);
                            },
                            _ => panic!("Invalid lookup data"),
                        },
                    }

                    match (seen, spec.allow_missing) {
                        (true, _) => (),
                        (false, true) => {
                            log::warn!("No data found for key {}", &buf);
                            continue;
                        },
                        (false, false) => {
                            return Err(format!("No data found for key {}", &id).into());
                        },
                    }

                    serde_json::to_writer(&mut writer, &data)?;
                    writeln!(&mut writer)?;
                    buf.clear();
                }
            },
        }
    }

    Ok(())
}
