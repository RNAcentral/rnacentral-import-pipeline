nextflow.enable.dsl=2

process slack_message {
    input:
    val(message)

    """
    #!/bin/bash
    curl -X POST -H 'Content-type: application/json' --data '{"text":"$message"}' $LITSCAN_SLACK_WEBHOOK
    """
}

process get_ids {
    input:
    path(database)

    output:
    tuple val("$database.baseName"), path('results')

    script:
    """
    psql -t -A -f $database "$PGDB_EMBASSY_USER" > results
    """
}

process check_ids {
    input:
    tuple val(database), file(results)

    output:
    tuple val(database), path('output')

    script:
    """
    litscan-check-ids.py $database $results output
    """
}

process sort_ids {
    publishDir "$baseDir/workflows/litscan/results/", mode: 'copy'

    input:
    tuple val(database), file(output)

    output:
    tuple val(database), path("${database}.txt")

    script:
    """
    LC_ALL=C sort -ufb $output > ${database}.txt
    """
}

process prepare_to_submit {
    memory '2GB'
    publishDir "$baseDir/workflows/litscan/submit/", mode: 'copy'

    input:
    tuple val(database), path("${database}.txt")

    output:
    tuple val(database), path("${database}_ids.txt")

    script:
    """
    litscan-get-unique-ids.sh ${database}.txt $database
    """
}

process submit_ids {
    input:
    val(_flag)

    output:
    val('done')

    script:
    """
    # create file with IDs (without URS)
    find $baseDir/workflows/litscan/submit/ -name "*.txt" -print0 | xargs -0 sed '/^URS/d' > all_ids.txt
    tr A-Z a-z < all_ids.txt | sort -u > all_ids_sorted.txt

    # get ids already scanned by LitScan (Expert DB ids only)
    psql -v ON_ERROR_STOP=1 -c "COPY (SELECT job_id FROM litscan_job WHERE job_id not like 'urs%' ORDER BY job_id) TO STDOUT;" $PGDATABASE > old_ids.txt
    sort old_ids.txt > old_ids_sorted.txt

    # get new ids
    comm -13 old_ids_sorted.txt all_ids_sorted.txt > new_ids.txt

    if [ -s new_ids.txt ]; then
      # new_ids.txt is not empty
      # get original ids (not in lowercase)
      # TODO: improve performance. This step may take a long time depending on the number of IDs.
      while IFS= read -r line; do
        grep -ixF "\$line" all_ids.txt | head -1 >> results.txt
      done < new_ids.txt

      # submit new ids only
      litscan-upload-ids.sh results.txt
      count=\$(sed -n '\$=' results.txt)
      curl -X POST -H 'Content-type: application/json' --data '{"text":"'\${count}' new id/gene/synonym submitted"}' \$LITSCAN_SLACK_WEBHOOK
    else
      # new_ids.txt is empty
      curl -X POST -H 'Content-type: application/json' --data '{"text":"No new id/gene/synonym to submit"}' $LITSCAN_SLACK_WEBHOOK
    fi
    """
}

process submit_urs {
    input:
    val(_flag)

    output:
    val('done')

    script:
    """
    # create file with URS only
    find $baseDir/workflows/litscan/submit/ -name "*.txt" -print0 | xargs -0 sed '/^URS/!d' > output.txt
    tr a-z A-Z < output.txt | sort -u > urs.txt

    # get URS already scanned by LitScan
    psql -v ON_ERROR_STOP=1 -c "COPY (SELECT job_id FROM litscan_job WHERE job_id like 'urs%' ORDER BY job_id) TO STDOUT;" $PGDATABASE > urs_lower.txt
    tr a-z A-Z < urs_lower.txt > urs_in_db.txt

    # get new URS
    comm -13 urs_in_db.txt urs.txt > new_urs.txt

    if [ -s new_urs.txt ]; then
      # new_urs.txt is not empty
      # submit new URS
      litscan-upload-ids.sh new_urs.txt
      count=\$(sed -n '\$=' new_urs.txt)
      curl -X POST -H 'Content-type: application/json' --data '{"text":"'\${count}' new URS submitted"}' \$LITSCAN_SLACK_WEBHOOK
    else
      # new_urs.txt is empty
      curl -X POST -H 'Content-type: application/json' --data '{"text":"No new URS to submit"}' $LITSCAN_SLACK_WEBHOOK
    fi
    """
}

workflow search_new_ids {
    emit: done
    main:
      Channel.of("Starting LitScan pipeline") | slack_message

      Channel.fromPath('workflows/litscan/queries/*.sql') \
      | get_ids \
      | check_ids \
      | sort_ids \
      | prepare_to_submit \
      | collect \
      | submit_ids \
      | submit_urs \
      | set { done }
}

workflow {
  search_new_ids()
}
