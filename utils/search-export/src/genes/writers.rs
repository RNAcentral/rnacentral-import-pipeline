use std::{
    collections::HashMap,
    fs::{
        create_dir_all,
        File,
    },
    io::{
        BufReader,
        BufWriter,
        Write,
    },
    path::{
        Path,
        PathBuf,
    },
};

use anyhow::{
    Context,
    Result,
};
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
            RegionGrouper,
            UrsRegion,
        },
    },
    search_xml::SearchEntry,
    sequences::{
        normalized::Normalized,
        so_tree,
    },
};

pub fn assembly_writer(base: &Path, assembly: &str) -> Result<BufWriter<File>> {
    let assembly_name = assembly.replace(".", "_");
    let mut path = PathBuf::from(base);
    path.push(assembly_name);
    path.set_extension("json");
    Ok(BufWriter::new(File::create(path)?))
}

pub fn write_split_selected(locus_path: &Path, sequence_file: &Path, output: &Path) -> Result<()> {
    create_dir_all(output).with_context(|| format!("Failed to create directory {output:?}"))?;

    let locus_file = File::open(locus_path)
        .with_context(|| format!("Failed to open locus file {locus_path:?}"))?;
    let locus_reader = BufReader::new(locus_file);
    let locus = Deserializer::from_reader(locus_reader).into_iter::<Grouped<UrsRegion>>();
    let locus = locus.map(|r| r.unwrap());
    let locus = locus.filter(|l| !l.is_empty());
    let locus = RegionGrouper::new(locus);
    let locus = locus.map(|l| l.unwrap());
    let locus = locus.map(|l| (*l.id(), l)).assume_sorted_by_key();

    let sequence_file = File::open(sequence_file)
        .with_context(|| format!("Failed to open sequence file {sequence_file:?}"))?;
    let sequence_reader = BufReader::new(sequence_file);
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
                .or_insert_with(|| assembly_writer(output, assembly_id).unwrap());
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

pub fn write_merged_members(so_file: &Path, member_file: &Path, output: &Path) -> Result<()> {
    let open_member = File::open(member_file)
        .with_context(|| format!("Failed to open member file {member_file:?}"))?;
    let member_reader = BufReader::new(open_member);
    let members = Deserializer::from_reader(member_reader).into_iter::<GeneMember>();
    let members = members.map(|r| r.unwrap());
    let mut members: Vec<GeneMember> = members.collect();
    members.sort_by_key(|m| m.gene_id());
    let grouped = members.into_iter().group_by(|l| l.gene_id());

    let tree = so_tree::load(so_file)
        .with_context(|| format!("Failed to parse SO tree file {so_file:?}"))?;

    let create =
        File::create(output).with_context(|| format!("Failed to create output file {output:?}"))?;
    let mut writer = BufWriter::new(create);
    for (locus_id, locus) in &grouped {
        log::info!("Merging locus id {locus_id}");
        let gene = Gene::new(locus, &tree)
            .with_context(|| format!("Failed to create gene for {locus_id}"))?;
        serde_json::to_writer(&mut writer, &gene)
            .with_context(|| format!("Failed to format to JSON: {gene:?}"))?;
        writeln!(&mut writer)?;
    }
    writer.flush()?;

    Ok(())
}

pub fn write_search_files(gene_file: &Path, xml_output: &Path, count_output: &Path) -> Result<()> {
    let open =
        File::open(gene_file).with_context(|| format!("Failed to open gene file {gene_file:?}"))?;

    let gene_reader = BufReader::new(open);
    let genes = Deserializer::from_reader(gene_reader).into_iter::<Gene>();

    let xml_file = File::create(xml_output)
        .with_context(|| format!("Failed to create XML output file {xml_output:?}"))?;
    let mut xml_writer = BufWriter::new(xml_file);
    let mut count = 0;
    let now: DateTime<Utc> = Utc::now();
    let date = now.format("%d/%m/%Y");
    writeln!(&mut xml_writer, "<database><name>RNAcentral</name><description></description><release>1.0</release><release_date>{}</release_date><entries>", &date)?;
    for gene in genes {
        let gene = gene.with_context(|| "Failed to deserialize a gene")?;
        let xml = gene.search_entry();
        let element = quick_xml::se::to_string(&xml)?;
        let element = element.replace(r"", "\\b");
        writeln!(&mut xml_writer, "{}", element)?;
        count += 1;
    }

    writeln!(&mut xml_writer, "</entries>")?;
    writeln!(&mut xml_writer, "<entry_count>{}</entry_count></database>", &count)?;
    xml_writer.flush()?;

    let mut count_writer = BufWriter::new(File::create(count_output)?);
    write!(&mut count_writer, "{}", &count)?;
    count_writer.flush()?;

    Ok(())
}
