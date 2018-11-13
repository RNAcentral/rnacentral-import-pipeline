#!/usr/bin/env nextflow

process species_to_map {
  executor 'local'

  input:
  file(query) from Channel.fromPath('files/genome-mapping/mappable.sql')

  output:
  stdout into raw_genomes
  // Writes: species, assembly, url_pattern
  // ex: homo_sapiens,CRGh38,url_pattern

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE"
  """
}

raw_genomes
  .splitCsv()
  .into { assemblies, genomes_to_fetch }

process fetch_unmapped_sequences {
  scratch true

  input:
  set val(species), val(assmebly), val(division) from assemblies

  output:
  set species, file('parts/*.fasta') into split_sequences

  script
  """
  set -o pipefail

  psql -v ON_ERROR_STOP=1 -v assembly_id=$assembly_id -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size ${params.qa.rfam_scan.chunk_size} --out-dir 'parts/' shuffled.fasta
  """
}

process download_genome {
  scratch true

  input:
  set val(species), val(assembly), val(division) from genomes_to_fetch

  output:
  set val(species), file('*.fa'), file('11.ooc') into genomes

  script:
  def url = GStringTemplateEngine()
    .createTemplate(params.genome_mapping.sources[division])
    .make([species: species])
    .toString()
  """
  fetch genome $url $species genome.fasta

  blat \
    -makeOoc=11.ooc \
    -stepSize=${params.genome_mapping.blat_options.step_size} \
    -repMatch=${params.genome_mapping.blat_options.rep_match} \
    -minScore=${params.genome_mapping.blat_options.min_score} \
    genome.fasta /dev/null /dev/null
  """
}

genomes
  .join(split_sequences)
  .flatMap { species, chrs, ooc_file, chunks ->
    [chrs, chunks].combinations().inject([]) { chr, chunk -> [species, ooc_file, chr, chunk] }
  }
  .set { targets }

process blat {
  input:
  set val(species), file(ooc), file(chromosome), file(chunk) from targets

  output:
  set file(species), file('output.psl') into blat_results

  """
  blat \
    -ooc=$ooc \
    -noHead \
    -q=rna \
    -stepSize=${params.genome_mapping.blat_options.step_size} \
    -repMatch=${params.genome_mapping.blat_options.rep_match} \
    -minScore=${params.genome_mapping.blat_options.min_score} \
    -minIdentity=${params.genome_mapping.blat_options.min_identity} \
    $chromosome $chunk output.psl
  """
}

 blat_results
  .groupTuple()
  .map { it -> it[1] }
  .set { species_results }

process select_mapped_locations {
  memory '10 GB'

  input:
  file('*.psl') from species_results

  output:
  file 'locations.csv' into selected_locations

  """
  set -o pipefail

  sort -t 10 *.psl | rnac genome-mapping select-hits - locations.csv
  """
}

process load_genome_mapping {
  maxForks 1

  input:
  file('raw*.csv') from selected_locations.collect()
  file(ctl) from Channel.fromPath('files/genome-mapping/load.ctl')

  output:
  val('done') into mapped_genomes

  script:
  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} genome-mapping
  """
}
