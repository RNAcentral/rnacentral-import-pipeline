#!/usr/bin/env nextflow

process fetch_crs_bed {
  input:
  val(remote) from Channel.from(params.crs.path)

  output:
  file('*.bed') into raw_crs mode flatten
  file("cfs.cfg") into raw_crs_config

  """
  cp "$remote/data/*.bed.gz" .
  gzip -d *.bed.gz
  cp "$remote/crs_filtered_bedfiles.cfg" "crs.cfg"
  """
}

raw_crs_config
  .splitCsv(sep: "\t")
  .map { row -> row[0] }
  .into { for_rfam_fetch; for_rnacentral_fetch }

raw_crs
  .map { f -> 
    def assembly =  file(f).name
      .replace("cmf_extend.", "")
      .replace(".fdr10.nonredundant.bed", "")
    [ assembly, f ]
  }
  .set { crs_bed_with_assemblies }

for_rfam_fetch
  .combine(Channel.fromPath("files/crs/fetch-rfam.sql"))
  .set { rfam_to_fetch }

process fetch_rfam_locations {
  input:
  set val(assembly), file(query) from rfam_to_fetch

  output:
  set val(assembly), file("${assembly}.bed") into rfam_coordinates

  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$assembly" -f $query "$PGDATABASE" > result.json
  rnac ftp-export coordinates as-bed result.json |\
  sort -k1,1 -k2,2n > ${assembly}.bed
  """
}

for_rnacentral_fetch
  .combine(Channel.fromPath("files/crs/fetch-rnacentral.sql"))
  .set { rnacentral_assemblies_to_fetch }

process fetch_rnacentral_bed {
  input:
  set val(assembly), file(query) from rnacentral_assemblies_to_fetch

  output:
  set val(assembly), file("{$assembly}.bed") into rnacentral_locations

  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$assembly" -f $query "$PGDATABASE" > result.json
  rnac ftp-export coordinates as-bed result.json |\
  sort -k1,1 -k2,2n > ${assembly}.bed
  """
}

crs_bed_with_assemblies
  .join(rfam_coordinates)
  .set { crs_to_clean }

process remove_rfam_crs {
  input:
  set val(assembly), file(crs), file(rfam) from crs_to_clean

  output:
  set val(assembly), file("cleaned.bed") into cleaned_crs

  """
  bedtools intersect -va $crs -b $rfam > cleaned.bed
  """
}

cleaned_crs
  .join(rnacentral_locations)
  .set { crs_to_select }

process select_rnacentral_crs {
  input:
  set val(assembly), file(crs), file(rnacentral_bed) from crs_to_select

  output:
  set val(assembly), file("selected.bed") into selected_crs

  """
  bedtools -a $crs -b $rnacentral > selected.bed
  """
}

process process_crs {
  input:
  set val(assembly), file(crs) from selected_crs

  output:
  file('complete_features.csv') into processed_crs

  """
  rnac crs $crs complete_features.csv
  """
}

process import_crs {
  input:
  file('complete_features*.csv') from processed_crs.collect()
  file(ctl) from Channel.fromPath('files/crs/load.ctl')

  """
  cp $ctl crs.ctl
  pgloader --on-error-stop crs.ctl
  """
}
