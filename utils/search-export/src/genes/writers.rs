use std::{
    collections::HashMap,
    fs::{
        create_dir_all,
        File,
    },
    io::{
        BufWriter,
        Write,
    },
    iter::FromIterator,
    path::{
        Path,
        PathBuf,
    },
};

use anyhow::Result;
use itertools::Itertools;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use rnc_core::{
    grouper::Grouped,
    psql::JsonlIterator,
};

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
    let mut path = PathBuf::from(base);
    path.push(assembly);
    path.set_extension("json");
    Ok(BufWriter::new(File::create(path)?))
}

pub fn write_split_selected(
    locus_file: &Path,
    sequence_file: &Path,
    output: &Path,
) -> Result<()> {
    create_dir_all(output)?;

    let locus: JsonlIterator<File, Grouped<Region>> = JsonlIterator::from_path(locus_file)?;
    let locus = locus.filter(|l| !l.is_empty());
    let locus = RegionGrouper::new(locus);
    let locus = locus.map(|l| (*l.id(), l)).assume_sorted_by_key();

    let sequences: JsonlIterator<File, Normalized> = JsonlIterator::from_path(sequence_file)?;
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
        }
    }

    Ok(())
}

pub fn write_merged_members(member_file: &Path, output: &Path) -> Result<()> {
    let members = JsonlIterator::from_path(member_file)?;
    let mut members: Vec<GeneMember> = members.collect();
    members.sort_by_key(|m| m.locus_id());
    let grouped = members.into_iter().group_by(|l| l.locus_id());

    let mut writer = BufWriter::new(File::create(&output)?);
    for (_locus_id, locus) in &grouped {
        let gene = Gene::from_iter(locus);
        serde_json::to_writer(&mut writer, &gene)?;
        writeln!(&mut writer)?;
    }

    Ok(())
}

pub fn write_search_files(gene_file: &Path, xml_output: &Path, count_output: &Path) -> Result<()> {
    let genes: JsonlIterator<File, Gene> = JsonlIterator::from_path(gene_file)?;

    let mut xml_writer = BufWriter::new(File::create(&xml_output)?);
    let mut count = 0;
    for gene in genes {
        let xml = gene.as_search();
        quick_xml::se::to_writer(&mut xml_writer, &xml)?;
        writeln!(&mut xml_writer)?;
        count += 1;
    }
    xml_writer.flush()?;

    let mut count_writer = BufWriter::new(File::create(&count_output)?);
    write!(&mut count_writer, "{}", &count)?;
    count_writer.flush()?;

    Ok(())
}
