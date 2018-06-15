md5_query = Channel.fromPath('files/ftp-export/md5/query.sql')
md5_readme = Channel.fromPath('files/ftp-export/md5/readme.txt')
process md5 {
  publishDir "${params.ftp_export}/md5/", mode: 'copy'

  input:
  file query from md5_query
  file 'readme.txt' form md5_readme

  output:
  file "example.txt" into __md5_files
  file "md5.tsv.gz" into __md5_files
  file "readme.txt" into __md5_files

  """
  psql -f "$query" "$PGDATABASE" > md5.tsv
  head -10 md5.tsv > example.tsv
  gzip md5.tsv
  """
}

id_mapping_query = Channel.fromPath('files/ftp-export/id-mapping/query.sql')
id_mapping_readme = Channel.fromPath('files/ftp-export/id-mapping/readme.txt')
process id_mapping {
  publishDir "${params.ftp_export}/id_maping/", mode: 'copy'

  input:
  file query from id_mapping_query
  file 'readme.txt' from id_mapping_readme

  output:
  file 'id_mapping.tsv.gz' into id_mapping
  file "example.txt" into __id_mapping_files
  file "readme.txt" into __id_mapping_files

  """
  psql -f "$query" "$PGDATABASE" | rnac ftp-export id-mapping - id_mapping.tsv
  head id_mapping.tsv > example.tsv
  gzip id_mapping.tsv
  """
}

process database_id_mapping {
  publishDir "${params.ftp_export}/id_maping/database_mappings/", mode: 'copy'

  input:
  file 'id_mapping.tsv.gz' into id_mapping

  output:
  file '*.tsv' into __database_mappings

  """
  zcat id_mapping.tsv.gz | awk '{ print >> ($1 ".tsv") }'
  """
}

rfam_annotations_query = Channel.from('files/ftp-export/rfam-annotations.sql')
process rfam_annotations {
  publishDir "${params.ftp_export}/rfam/", mode: 'copy'

  input:
  file query from rfam_annotations_query

  output:
  file 'rfam_annotations.tsv.gz' into __rfam_files
  file "example.txt" into __rfam_files
  file "readme.txt" into __rfam_files

  """
  psql -f "$query" "$PGDATABASE" | gzip > rfam_annotations.tsv.gz
  """
}

process inactive_fasta {
  publishDir "${params.ftp_export}/sequences/", mode: 'copy'

  output:
  file 'rnacentral_inactive.fasta.gz' into __sequences

  """
  rnac ftp-export sequences inactive - | gzip > rnacentral_inactive.fasta.gz
  """
}

process active_fasta {
  publishDir "${params.ftp_export}/sequences/", mode: 'copy'

  output:
  file 'rnacentral_active.fasta.gz' into active_sequences
  file 'example.txt' into __sequences
  file 'readme.txt' into __sequences

  """
  rnac ftp-export sequences active rnacentral_active.fasta
  rnac ftp-export sequences readme readme.txt
  head rnacentral_active.fasta > example.txt
  gzip rnacentral_active.fasta
  """

}

process species_specific_fasta {
  publishDir "${params.ftp_export}/sequences/", mode: 'copy'

  output:
  file 'rnacentral_species_specific_ids.fasta.gz' into __sequences

  """
  rnac ftp-export sequences species-specific - | rnacentral_species_specific_ids.fasta.gz
  """
}

process nhmmer_fasta {
  publishDir "${params.ftp_export}/sequences/.internal/", mode: 'copy'

  input:
  file 'rnacentral_active.fasta.gz' from active_sequences

  output:
  file 'rnacentral_nhmmer_*.fasta' into __sequences

  """
  zcat rnacentral_active.fasta.gz | rnac ftp-export sequences split-nhmmer -
  """
}

process find_ensembl_chunks {
  output:
  stdout raw_ensembl_ranges

  """
  rnac ftp-export ensembl ranges
  """
}

ensembl_ranges = raw_ensembl_ranges.splitCsv(sep='\t')

process ensembl_export_chunk {
  publishDir "${params.ftp_export}/json/", mode: 'copy'

  input:
  set val(min), val(max) from ensembl_ranges

  output:
  file "ensembl-xref-$min-${max}.json" into __ensembl_export

  """
  rnac ftp-export ensembl $min $max ensembl-xref-$min-${max}.json
  """
}

process go_annotation {
  publishDir "${params.ftp_export}/go_annotations/", mode: 'copy'

  input:
  file query from go_annotation_query

  output:
  file "rnacentral_rfam_annotations.tsv.gz" into __go_annotations

  """
  psql -f "$query" "$PGDATABASE" | rnac ftp-export rfam-go-annotations - - | gzip > rnacentral_rfam_annotations.tsv.gz
  """
}

process genome_coordiantes {
  publishDir "${params.ftp_export}/genome_coordinates/", mode: 'copy'

  output:
  stdout raw_coordinates

  """
  rnac ftp-export has-coordinates
  """
}

raw_coordinates.splitCsv(sep='\t').into { coordinates }

process format_coordinates {
  publishDir "${params.ftp_export}/genome_coordinates/", mode: 'copy'

  input:
  set val format, val assembly, val species from coordinates

  output:
  set val(format), file(result) into formatted_coordinates

  script:
  result = "${species}.${assembly}.${format}.gz"
  """
  rnac ftp-export coordinates $format $species $assembly - | gzip > $result
  """
}


// formatted_coordiantes
//   .filter(it[0] == 'bed')
//   .map { it[1] }
//   .into { bed_files }

// process big_bed {
//   input:
//   file(bed_file) from bed_files

//   output:
//   file(result) into __big_bed

//   script:
//   result = "${bed_file.baseName}.bigBed"
//   """
//   wget -o chrom.size http://hgdownload.soe.ucsc.edu/goldenPath/$ucsc_assembly/bigZips/${ucsc_assembly}.chrom.sizes
//   sort -k1,1 -k2,2n $bed_file > sorted.bed
//   bedToBigBed sorted.bed chrom.size $result
//   """
// }

// process track_hub {
//   publishDir "${params.ftp_export}/genome_coordinates/track_hub/", mode: 'copy'

//   output:
//   stdout raw_trackhub_assemblies

//   """
//   rnac ftp-export has-trackhub
//   """
// }

// track_hub_assemblies = raw_trackhub_assemblies.splitCsv(sep='\t')

// process assembly_track_hub {
//   publishDir "${params.ftp_export}/genome_coordinates/track_hub/$assembly/", mode: 'copy'

//   input:
//   val assembly from track_hub_assemblies

//   """
//   rnac ftp-export track-hub $assembly
//   """
// }

// process gpi {
//   publishDir "${params.ftp_export}/gpi/", mode: 'copy'

//   output:
//   file 'rnacentral.gpi.gz' into __gpi

//   """
//   rnac ftp-export gpi - | gzip > rnacentral.gpi.gz
//   """
// }
