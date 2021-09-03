use std::{collections::HashMap, fs::{
        create_dir_all,
        File,
    }, io::{BufReader, BufWriter, Write}, iter::FromIterator, path::{
        Path,
        PathBuf,
    }};

use anyhow::Result;
use chrono::{
    DateTime,
    Utc,
};
use itertools::Itertools;

use serde_json::Deserializer;
use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use rnc_core::grouper::Grouped;

use crate::{
    genes::{
        gene::Gene,
        gene_member::GeneMember,
        region::{
            Region,
            RegionGrouper,
        },
    },
    sequences::normalized::Normalized,
};

pub fn assembly_writer(base: &Path, assembly: &str) -> Result<BufWriter<File>> {
    let assembly_name = assembly.replace(".", "_");
    let mut path = PathBuf::from(base);
    path.push(assembly_name);
    path.set_extension("json");
    Ok(BufWriter::new(File::create(path)?))
}

pub fn write_split_selected(locus_path: &Path, sequence_file: &Path, output: &Path) -> Result<()> {
    create_dir_all(output)?;

    let locus_reader = BufReader::new(File::open(locus_path)?);
    let locus = Deserializer::from_reader(locus_reader).into_iter::<Grouped<Region>>();
    let locus = locus.map(|r| r.unwrap());
    let locus = locus.filter(|l| !l.is_empty());
    let locus = RegionGrouper::new(locus);
    let locus = locus.map(|l| (*l.id(), l)).assume_sorted_by_key();

    let sequence_reader = BufReader::new(File::open(sequence_file)?);
    let sequences = Deserializer::from_reader(sequence_reader).into_iter::<Normalized>();
    let sequences = sequences.map(|r| r.unwrap());
    let sequences = sequences.map(|n| (*n.id(), n)).assume_sorted_by_key();

    let mut writers = HashMap::new();
    let merged = locus.join(sequences);

    for (_id, pair) in merged {
        let (seq_locus, sequence) = pair;
        for (assembly_id, regions) in seq_locus.regions_by_assembly() {
            let mut writer = writers
                .entry(assembly_id.clone())
                .or_insert_with(|| assembly_writer(&output, &assembly_id).unwrap());
            for region in regions {
                let member = GeneMember::new(region.clone(), sequence.clone());
                serde_json::to_writer(&mut writer, &member)?;
                writeln!(&mut writer)?;
            }
            writer.flush()?;
        }
    }

    for (_assembly, writer) in writers.iter_mut() {
        writer.flush()?;
    }

    Ok(())
}

pub fn write_merged_members(member_file: &Path, output: &Path) -> Result<()> {
    let member_reader = BufReader::new(File::open(member_file)?);
    let members = Deserializer::from_reader(member_reader).into_iter::<GeneMember>();
    let members = members.map(|r| r.unwrap());
    let mut members: Vec<GeneMember> = members.collect();
    members.sort_by_key(|m| m.locus_id());
    let grouped = members.into_iter().group_by(|l| l.locus_id());

    let mut writer = BufWriter::new(File::create(&output)?);
    for (_locus_id, locus) in &grouped {
        let gene = Gene::from_iter(locus);
        serde_json::to_writer(&mut writer, &gene)?;
        writeln!(&mut writer)?;
    }
    writer.flush()?;

    Ok(())
}

pub fn write_search_files(gene_file: &Path, xml_output: &Path, count_output: &Path) -> Result<()> {

    let gene_reader = BufReader::new(File::open(gene_file)?);
    let genes = Deserializer::from_reader(gene_reader).into_iter::<Gene>();

    let mut xml_writer = BufWriter::new(File::create(&xml_output)?);
    let mut count = 0;
    let now: DateTime<Utc> = Utc::now();
    let date = now.format("%d/%m/%Y");
    writeln!(&mut xml_writer, "<database><name>RNAcentral</name><description></description><release>1.0</release><release_date>{}</release_date><entries>", &date)?;
    for gene in genes {
        let mut gene = gene?;
        gene.fill_default_name();
        let xml = gene.as_search();
        let element = quick_xml::se::to_string(&xml)?;
        let element = element.replace(r"", "\\b");
        writeln!(&mut xml_writer, "{}", element)?;
        count += 1;
    }

    writeln!(&mut xml_writer, "</entries>")?;
    writeln!(&mut xml_writer, "<entry_count>{}</entry_count></database>", &count)?;
    xml_writer.flush()?;

    let mut count_writer = BufWriter::new(File::create(&count_output)?);
    write!(&mut count_writer, "{}", &count)?;
    count_writer.flush()?;

    Ok(())
}
