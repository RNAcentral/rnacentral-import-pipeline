if [ -z "$1" ]
  then
    echo "Specify location of psl files"
    return
fi

if [ -z "$2" ]
  then
    echo "Specify assembly id"
    return
fi

INPUT_FOLDER=$1
ASSEMBLY_ID=$2
echo 'Processing files in ' $1

MATCHES_PSL=${INPUT_FOLDER}/matches.psl
MATCHES_BED=${INPUT_FOLDER}/matches-split.bed

for i in {1..2}; do
	if [ $i -eq 1 ]; then
	    OUTPUT=${INPUT_FOLDER}/mapping.tsv
	else
	    OUTPUT=${INPUT_FOLDER}/mapping-inexact.tsv
	fi

	echo '' > $MATCHES_PSL
	for d in $INPUT_FOLDER/*.psl; do

		  TEMPFILE=tempfile.psl
		    if [ $i -eq 1 ]; then
				awk '{
				  matches = $1
				  seqlength = $11
				  if (matches == seqlength) { print }
				}' $d > $TEMPFILE
		    elif [ $i -eq 2 ]; then
		        awk '{
		          if ($11>15 && $1/$11>0.95 && $1/$11<1) { print }
		        }' $d > $TEMPFILE
		    fi

			cat $TEMPFILE | \
			# skip long introns in short sequences
			awk '{
			  matches = $1
			  target_insertions = $8
			  if (matches < 100 && target_insertions > 10) {;} else  { print }
			}' | \
			# print out modified psl with updated strand and region_id
			awk -v OFS='\t' '{
			  psl_strand = $9
			  if (psl_strand == "+") strand=1;
			  else if (psl_strand == "-") strand=-1;
			  else strand = 0;
			  start = $16+1;
			  region_id = $10 "@" $14 "/" start "-" $17 ":" $9
			  print $1, $2, $3, $4, $5, $6, $7, $8, strand, region_id, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21
			}' >> $MATCHES_PSL
	done

	# use bedops to convert psl to bed, use --split to convert bed blocks into separate lines for each exon
	psl2bed --split < $MATCHES_PSL > $MATCHES_BED

	# transform bed to tsv that can be loaded in the database
	awk -v OFS='\t' '{
	  identity = $5/$7
	  print $4, $1, $2, $3, $6, identity
	}' $MATCHES_BED | \
	awk -v OFS='\t' -v assembly_id="$ASSEMBLY_ID" '{
	  split($1, a, "@");
	  split(a[1], b, "_");
	  start = $3+1
	  print a[1], b[1], b[2], c[1], $1, $2, start, $4, $5, assembly_id, $6, NR
	}' > $OUTPUT

	rm $MATCHES_PSL
	rm $MATCHES_BED
  rm $TEMPFILE

	echo 'See file ' $OUTPUT

done
