process find_genome_coordinate_jobs {
  input:
  file query from Channel.fromPath('files/ftp-export/genome_coordinates/known-coordinates.sql')

  output:
  file('coordinates.txt') into species_to_format

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > coordinates.txt
  """
}

species_to_format
  .splitCsv()
  .combine(Channel.fromPath('files/ftp-export/genome_coordinates/query.sql'))
  .set { coordinates_to_fetch }

process coordinate_readme {
  publishDir "${params.ftp_export.publish}/genome_coordinates/", mode: 'copy'

  input:
  file raw from Channel.fromPath('files/ftp-export/genome_coordinates/readme.mkd')

  output:
  file('readme.txt') into __coordinate_readmes

  """
  cat $raw > readme.txt
  """
}

process fetch_raw_coordinate_data {
  maxForks params.ftp_export.coordinate.maxForks

  input:
  set val(assembly), val(species), val(taxid), file(query) from coordinates_to_fetch

  output:
  set val(assembly), val(species), file('result.json') into raw_coordinates

  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$assembly" -f $query "$PGDATABASE" > result.json
  """
}

raw_coordinates.into { bed_coordinates; gff_coordinates }

process format_bed_coordinates {
  publishDir "${params.ftp_export.publish}/genome_coordinates/bed/", mode: 'copy'

  input:
  set val(assembly), val(species), file(raw_data) from bed_coordinates

  output:
  set val(assembly), file(result) into bed_files

  script:
  result = "${species}.${assembly}.bed.gz"
  """
  set -o pipefail

  rnac ftp-export coordinates as-bed $raw_data |\
  sort -k1,1 -k2,2n |\
  gzip > $result
  """
}

// process generate_big_bed {
//   publishDir "${params.ftp_export.publish}/genome_coordinates/bed", mode: 'copy'

//   input:
//   set val assembly, file bed_file from bed_files

//   output:
//   file(bigBed) into big_bed_files

//   script:
//   chrom = "${bed_file.baseName}.chrom.sizes"
//   bigBed = "${bed_file.baseName}.bigBed"
//   """
//   fetchChromSizes "$assembly" > "$chrom"
//   bedToBigBed -type -type=bed12+3 bed_file $chrom > $bigBed
//   """
// }

process generate_gff3 {
  memory '8 GB'

  publishDir "${params.ftp_export.publish}/genome_coordinates/gff3", mode: 'move'

  input:
  set val(assembly), val(species), file(raw_data) from gff_coordinates

  output:
  file result into gff3_files

  script:
  result = "${species}.${assembly}.gff3.gz"
  """
  rnac ftp-export coordinates as-gff3 $raw_data - | gzip > $result
  """
}
