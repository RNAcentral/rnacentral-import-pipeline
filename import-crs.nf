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
  set val(assembly), file("selected_crs.tsv") into selected_crs

  shell:
  '''
  bedtools intersect -a !{crs} -b !{rfam} -v > crs-no-rfam.bed
  bedtools sort -i crs-no-rfam.bed > cleaned.bed

  bed12ToBed6 -i !{rnacentral} |\
  awk -n 'BEGIN { OFS="\t" } { print $1,$2,$3,$4":"$1":"$2":"$3":"$11,$5,$6,$7,$8,$9,$10,$11,$12 }' |\
  awk 'BEGIN { OFS="\t" } {
    split($4,a,":");
    split(a[length(a)],b,",");
    offset=0;
    if($5!=1) {
      for(i=1;i<$5;i++) {
        if ($6=="+") {
          offset+=b[i]
        } else {
          offset=offset+b[length(b)-i+1]
        }
      }
    };
    print $1,$2,$3,$4":"offset,$5,$6
  }' > unsorted-exons.bed

  bedtools sort -i unsorted-exons.bed > exons.bed

  bedtools intersect -a cleaned.bed -b exons.bed -wo > intersect_exon.bed
  awk 'BEGIN {
    OFS="\t";
    print "URS_taxid","CRS_start_relative_to_URS","CRS_end_relative_to_URS","chromosome","CRS_start_relative_to_genome","CRS_end_relative_to_genome","CRS_id","CRS_fdr"
  } {
    split($10,a,":");
  if ($12=="+") {
    start=$2-$8+1+a[$6];
    end=$3-$8+a[$6]
  }
  else {
    start=$9-$3+1+a[$6];
    end=$9-$2+a[$6]
  };
  print a[1],start,end,$1,$2+1,$3,$4,$5
  }' intersect_exon.bed > selected_crs.tsv
  '''
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
