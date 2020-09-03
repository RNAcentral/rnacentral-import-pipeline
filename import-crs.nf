#!/usr/bin/env nextflow

process fetch_crs_bed {
  input:
  val(remote) from Channel.from(params.crs.path)

  output:
  file('*.bed') into raw_crs mode flatten

  """
  cp $remote/data/*.bed.gz .
  gzip -d *.bed.gz
  """
}

raw_crs
  .map { f ->
    def assembly =  file(f).name
      .replace("cmf_extend.", "")
      .replace(".fdr10.nonredundant.bed", "")
    [ assembly, f ]
  }
  .into { crs_bed_with_assemblies; pre_fetch }

pre_fetch
  .map { assembly, _ -> [assembly, params.crs.assembly_rnacentral_mapping[assembly]] }
  .into { for_rfam_fetch; for_rnacentral_fetch }

for_rfam_fetch
  .combine(Channel.fromPath("files/crs/fetch-rfam.sql"))
  .set { rfam_to_fetch }

process fetch_rfam_locations {
  maxForks 4

  input:
  set val(crs_assembly), val(assembly), file(query) from rfam_to_fetch

  output:
  set val(crs_assembly), file("rfam-${assembly}.bed") into rfam_coordinates

  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$crs_assembly" -f $query "$PGDATABASE" > result.json
  rnac ftp-export coordinates as-bed result.json result.bed
  bedtools sort -i result.bed > rfam-${assembly}.bed
  """
}

for_rnacentral_fetch
  .combine(Channel.fromPath("files/crs/fetch-rnacentral.sql"))
  .set { rnacentral_assemblies_to_fetch }

process fetch_rnacentral_bed {
  maxForks 4

  input:
  set val(crs_assembly), val(assembly), file(query) from rnacentral_assemblies_to_fetch

  output:
  set val(crs_assembly), file("rnacentral-${assembly}.bed") into rnacentral_locations

  """
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$assembly" -f $query "$PGDATABASE" > result.json
  rnac ftp-export coordinates as-bed result.json > result.bed
  bedtools sort -i result.bed > rnacentral-${assembly}.bed
  """
}

crs_bed_with_assemblies
  .join(rfam_coordinates)
  .join(rnacentral_locations)
  .set { crs_to_clean }

process remove_rfam_crs {
  input:
  set val(assembly), file(crs), file(rfam), file(rnacentral) from crs_to_clean

  output:
  file('complete_features.csv') into processed_crs

  script:
  def must_clean = params.crs.must_clean_bed.contains(assembly) ? '1' : '0';
  """
  crs-overlaps $crs $rfam $rnacentral $must_clean
  rnac crs selected_crs.tsv complete_features.csv
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
