#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process setup {
  when { params.genome_mapping.run }

  input:
  val(_flag)
  path(query)

  output:
  path('species.csv')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > species.csv
  """
}

process fetch_unmapped_sequences {
  tag { species }
  maxForks params.genome_mapping.fetch_unmapped_sequences.directives.maxForks
  memory '10GB'
  clusterOptions '-sp 100'

  input:
  tuple val(species), val(assembly_id), val(taxid), val(division), path(possible_query), path(mapped_query), path(query)

  output:
  tuple val(species), path('parts/*.fasta')

  """
  set -euo pipefail

  psql \
    -v ON_ERROR_STOP=1 \
    -v taxid=$taxid \
    -v min_length=${params.genome_mapping.min_length} \
    -v max_length=${params.genome_mapping.max_length} \
    -f "$possible_query" \
    "$PGDATABASE" | sort > possible

  psql \
    -v ON_ERROR_STOP=1 \
    -v assembly_id=$assembly_id \
    -f "$mapped_query" \
    "$PGDATABASE" | sort -u > mapped

  mkdir parts

  comm -23 possible mapped | awk 'BEGIN { FS="_"; OFS="," } { print \$1, \$0 }' > urs-to-compute
  if [[ -s urs-to-compute ]]; then
    psql -q -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
    json2fasta raw.json rnacentral.fasta
    seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta

    split-sequences \
      --max-nucleotides ${params.genome_mapping.fetch_unmapped_sequences.nucleotides_per_chunk} \
      --max-sequences ${params.genome_mapping.fetch_unmapped_sequences.sequences_per_chunk} \
        shuffled.fasta parts
  else
    touch parts/empty.fasta
  fi
  """
}

process download_genome {
  tag { species }
  memory { params.genome_mapping.download_genome.directives.memory }
  errorStrategy 'ignore'

  input:
  tuple val(species), val(assembly), val(taxid), val(division)

  output:
  tuple val(species), val(assembly), path("${species}.{2bit,ooc}")

  """
  set -o pipefail

  rnac genome-mapping url-for --host=$division $species $assembly - |\
    xargs -I {} fetch generic '{}' ${species}.fasta.gz

  gzip -d ${species}.fasta.gz
  faToTwoBit -noMask ${species}.fasta ${species}.2bit
  blat \
    -makeOoc=${species}.ooc \
    -stepSize=${params.genome_mapping.blat.options.step_size} \
    -repMatch=${params.genome_mapping.blat.options.rep_match} \
    -minScore=${params.genome_mapping.blat.options.min_score} \
    ${species}.fasta /dev/null /dev/null
  """
}

process blat {
  tag { "${species}-${genome.baseName}-${chunk.baseName}" }
  memory { params.genome_mapping.blat.directives.memory }
  errorStrategy 'ignore'

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

  output:
  val('done')

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} genome-mapping
  split-and-load $attempted_ctl 'attempted*.csv' ${params.import_data.chunk_size} genome-mapping-attempted
  """
}

workflow genome_mapping {
  take: ready
  emit: done
  main:
    Channel.fromPath('files/genome-mapping/find_species.sql').set { find_species }
    Channel.fromPath('files/genome-mapping/possible.sql').set { possible_sql }
    Channel.fromPath('files/genome-mapping/get-mapped.sql').set { mapped_sql }
    Channel.fromPath('files/genome-mapping/find-unmapped.sql').set { unmapped_sql }
    Channel.fromPath('files/genome-mapping/load.ctl').set { hits_ctl }
    Channel.fromPath('files/genome-mapping/attempted.ctl').set { attempted_ctl }

    setup(ready, find_species) \
    | splitCsv \
    | filter { s, a, t, d -> !params.genome_mapping.species_excluded_from_mapping.contains(s) } \
    | set { genome_info }

    genome_info \
    | combine(possible_sql) \
    | combine(mapped_sql) \
    | combine(unmapped_sql) \
    | fetch_unmapped_sequences \
    | set { split_sequences }

    genome_info \
    | download_genome \
    | join(split_sequences) \
    | flatMap { species, assembly, genome_chunks, chunks ->
      [genome_chunks.collate(2), chunks]
        .combinations()
        .inject([]) { acc, files -> acc << [species, assembly] + files.flatten() }
    } \
    | filter { s, a, g, o, chunk -> !chunk.isEmpty() } \
    | blat

    blat.out.hits | groupTuple | select_mapped_locations | collect | set { hits }
    blat.out.attempted | collect | set { attempted }

    load_mapping(hits, hits_ctl, attempted, attempted_ctl) | set { done }
}

workflow {
  genome_mapping(Channel.from('ready'))
}
