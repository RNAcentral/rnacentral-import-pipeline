#!/usr/bin/env nextflow

process fetch_crs_bed {
  input:
  val(remote)

  output:
  path('*.bed')

  script:
  """
  cp $remote/data/*.bed.gz .
  gzip -d *.bed.gz
  """
}

process fetch_rfam_locations {
  maxForks 4

  input:
  tuple val(crs_assembly), val(assembly), path(query)

  output:
  tuple val(crs_assembly), path("rfam-${assembly}.bed")

  script:
  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$crs_assembly" -f $query "\$PGDATABASE" > result.json
  rnac ftp-export coordinates as-bed result.json result.bed
  bedtools sort -i result.bed > rfam-${assembly}.bed
  """
}

process fetch_rnacentral_bed {
  maxForks 4

  input:
  tuple val(crs_assembly), val(assembly), path(query)

  output:
  tuple val(crs_assembly), path("rnacentral-${assembly}.bed")

  script:
  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$assembly" -f $query "\$PGDATABASE" > result.json
  rnac ftp-export coordinates as-bed result.json > result.bed
  bedtools sort -i result.bed > rnacentral-${assembly}.bed
  """
}

process find_rnacentral_crs_features {
  input:
  tuple val(assembly), path(crs), path(rfam), path(rnacentral)

  output:
  path('complete_features.csv')

  script:
  def must_clean = params.crs.must_clean_bed.contains(assembly) ? '1' : '0';
  """
  crs-overlaps $crs $rfam $rnacentral $must_clean
  rnac crs selected_crs.tsv complete_features.csv
  """
}

process import_crs {
  input:
  path('complete_features*.csv')
  path(ctl)

  script:
  """
  cp $ctl crs.ctl
  pgloader --on-error-stop crs.ctl
  """
}

workflow {
  crs_bed_with_assemblies = fetch_crs_bed(channel.of(params.crs.path))
    | flatten
    | map { f ->
        def assembly = file(f).name
          .replace("cmf_extend.", "")
          .replace(".fdr10.nonredundant.bed", "")
        [ assembly, f ]
      }

  pre_fetch = crs_bed_with_assemblies
    .map { assembly, _f -> [assembly, params.crs.assembly_rnacentral_mapping[assembly]] }

  rfam_coordinates = fetch_rfam_locations(
    pre_fetch.combine(channel.fromPath("files/crs/fetch-rfam.sql"))
  )

  rnacentral_locations = fetch_rnacentral_bed(
    pre_fetch.combine(channel.fromPath("files/crs/fetch-rnacentral.sql"))
  )

  crs_to_clean = crs_bed_with_assemblies
    .join(rfam_coordinates)
    .join(rnacentral_locations)

  processed_crs = find_rnacentral_crs_features(crs_to_clean)

  import_crs(processed_crs.collect(), channel.fromPath('files/crs/load.ctl'))
}
