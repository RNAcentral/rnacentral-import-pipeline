#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process setup {
  input:
  path(setup)
  path(query)

  output:
  path('species.csv')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=${params.genome_mapping.to_map_table} \
    -v species_to_map=${params.genome_mapping.species_table} \
    -v min_length=${params.genome_mapping.min_length} \
    -v max_length=${params.genome_mapping.max_length} \
    -f "$setup" "$PGDATABASE"
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=${params.genome_mapping.to_map_table} \
    -v species_to_map=${params.genome_mapping.species_table} \
    -f "$query" "$PGDATABASE" > species.csv
  """
}

process fetch_unmapped_sequences {
  tag { species }
  maxForks params.genome_mapping.fetch_unmapped_sequences.directives.maxForks
  clusterOptions '-sp 100'

  input:
  tuple val(species), val(assembly_id), val(taxid), val(division), path(query)

  output:
  tuple val(species), path('parts/*.fasta')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v taxid=$taxid \
    -v assembly_id=$assembly_id \
    -v tablename=${params.genome_mapping.to_map_table} \
    -v species_to_map=${params.genome_mapping.species_table} \
    -f "$query" \
    "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta

  split-sequences \
    --max-nucleotides ${params.genome_mapping.fetch_unmapped_sequences.nucleotides_per_chunk} \
    --max-sequences ${params.genome_mapping.fetch_unmapped_sequences.sequences_per_chunk} \
      shuffled.fasta parts
  """
}

process download_genome {
  tag { species }
  memory { params.genome_mapping.download_genome.directives.memory }
  errorStrategy 'ignore'

  input:
  tuple val(species), val(assembly), val(taxid), val(division)

  output:
  tuple val(species), val(assembly), path('parts/*.{2bit,ooc}')

  """
  set -o pipefail

  rnac genome-mapping url-for --host=$division $species $assembly - |\
    xargs -I {} fetch generic '{}' ${species}.fasta.gz

  gzip -d ${species}.fasta.gz
  split-sequences \
    --max-nucleotides ${params.genome_mapping.download_genome.nucleotides_per_chunk} \
    --max-sequences ${params.genome_mapping.download_genome.sequences_per_chunk} \
    ${species}.fasta parts

  find parts -name '*.fasta' |\
    xargs -I {} faToTwoBit -noMask {} {}.2bit

  find parts -name '*.fasta' |\
  xargs -I {} \
    blat \
      -makeOoc={}.ooc \
      -stepSize=${params.genome_mapping.blat.options.step_size} \
      -repMatch=${params.genome_mapping.blat.options.rep_match} \
      -minScore=${params.genome_mapping.blat.options.min_score} \
      {} /dev/null /dev/null
  """
}

process blat {
  tag { "${species}-${genome.baseName}-${chunk.baseName}" }
  memory { params.genome_mapping.blat.directives.memory }
  errorStrategy 'finish'

  input:
  tuple val(species), val(assembly), path(genome), path(ooc), path(chunk)

  output:
  tuple val(species), path('selected.json'), emit: hits
  path 'attempted.csv', emit: attempted

  """
  set -o pipefail

  blat \
    -ooc=$ooc \
    -noHead \
    -q=rna \
    -stepSize=${params.genome_mapping.blat.options.step_size} \
    -repMatch=${params.genome_mapping.blat.options.rep_match} \
    -minScore=${params.genome_mapping.blat.options.min_score} \
    -minIdentity=${params.genome_mapping.blat.options.min_identity} \
    $genome $chunk output.psl

  sort -k 10 output.psl |\
    rnac genome-mapping blat serialize $assembly - - |\
    rnac genome-mapping blat select - selected.json

  rnac genome-mapping create-attempted $chunk $assembly attempted.csv
  """
}

process select_mapped_locations {
  tag { species }
  memory { '20GB' }

  input:
  tuple val(species), path('selected*.json')

  output:
  path('locations.csv')

  """
  set -o pipefail

  find . -name 'selected*.json' |\
    xargs cat |\
    rnac genome-mapping blat select --sort - - |\
    rnac genome-mapping blat as-importable - locations.csv
  """
}

process load_mapping {
  maxForks 1

  input:
  path('raw*.csv')
  path(ctl)
  path('attempted*.csv')
  path(attempted_ctl)

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} genome-mapping
  split-and-load $attempted_ctl 'attempted*.csv' ${params.import_data.chunk_size} genome-mapping-attempted
  psql -v ON_ERROR_STOP=1 -c 'DROP TABLE ${params.genome_mapping.to_map_table}' $PGDATABASE
  """
}

workflow genome_mapping {
  Channel.fromPath('files/genome-mapping/setup.sql').set { setup_sql }
  Channel.fromPath('files/genome-mapping/find-species.sql').set { find_species }
  Channel.fromPath('files/genome-mapping/load.ctl').set { hits_ctl }
  Channel.fromPath('files/genome-mapping/attempted.ctl').set { attempted_ctl }

  setup(setup_sql, find_species) \
  | splitCsv \
  | filter { s, a, t, d -> !params.genome_mapping.species_excluded_from_mapping.contains(s) } \
  | set { genome_info }

  genome_info \
  | download_genome \
  | set { genomes }

  genome_info \
  | combine(Channel.fromPath('files/genome-mapping/find-unmapped.sql')) \
  | fetch_unmapped_sequences \
  | set { split_sequences }

  genomes \
  | join(split_sequences) \
  | flatMap { species, assembly, genome_chunks, chunks ->
    [genome_chunks.collate(2), chunks]
      .combinations
      .inject([]) { acc, files -> acc << [species, assembly] + files.flatten() }
  } \
  | blat

  blat.out.hits \
  | groupTuple \
  | select_mapped_locations \
  | collect \
  | set { hits }

  blat.out.attempted | collect | set { attempted }

  load_mapping(hits, hits_ctl, attempted, attempted_ctl)
}

workflow {
  genome_mapping()
}
