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

process register_ids {
    input:
    val(_flag)

    output:
    path("registered_ids.txt")

    script:
    """
    # create file with IDs (without URS)
    find $baseDir/workflows/litscan/submit/ -name "*.txt" -print0 | xargs -0 sed '/^URS/d' > all_ids.txt
    tr A-Z a-z < all_ids.txt | sort -u > all_ids_sorted.txt

    # get ids already registered in LitScan (Expert DB ids only)
    psql -v ON_ERROR_STOP=1 -c "COPY (SELECT job_id FROM litscan_job WHERE job_id not like 'urs%' ORDER BY job_id) TO STDOUT;" \$PGDATABASE > old_ids.txt
    sort old_ids.txt > old_ids_sorted.txt

    # get new ids
    comm -13 old_ids_sorted.txt all_ids_sorted.txt > new_ids.txt

    touch registered_ids.txt
    if [ -s new_ids.txt ]; then
      # get original ids (not in lowercase) - single pass, no subprocess per line
      awk '
        NR==FNR { ids[tolower(\$0)] = 1; next }
        tolower(\$0) in ids && !seen[tolower(\$0)]++ { print }
      ' new_ids.txt all_ids.txt >> results.txt

      # register new ids in the database
      litscan-register-ids.py results.txt registered_ids.txt
      count=\$(wc -l < registered_ids.txt)
      curl -X POST -H 'Content-type: application/json' --data '{"text":"'\${count}' new id/gene/synonym registered for scanning"}' \$LITSCAN_SLACK_WEBHOOK
    else
      curl -X POST -H 'Content-type: application/json' --data '{"text":"No new id/gene/synonym to register"}' \$LITSCAN_SLACK_WEBHOOK
    fi
    """
}

process register_urs {
    input:
    path(registered_ids)

    output:
    path("registered_urs.txt")

    script:
    """
    # create file with URS only
    find $baseDir/workflows/litscan/submit/ -name "*.txt" -print0 | xargs -0 sed '/^URS/!d' > output.txt
    tr a-z A-Z < output.txt | sort -u > urs.txt

    # get URS already registered in LitScan
    psql -v ON_ERROR_STOP=1 -c "COPY (SELECT job_id FROM litscan_job WHERE job_id like 'urs%' ORDER BY job_id) TO STDOUT;" \$PGDATABASE > urs_lower.txt
    tr a-z A-Z < urs_lower.txt > urs_in_db.txt

    # get new URS
    comm -13 urs_in_db.txt urs.txt > new_urs.txt

    touch registered_urs.txt
    if [ -s new_urs.txt ]; then
      litscan-register-ids.py new_urs.txt registered_urs.txt
      count=\$(wc -l < registered_urs.txt)
      curl -X POST -H 'Content-type: application/json' --data '{"text":"'\${count}' new URS registered for scanning"}' \$LITSCAN_SLACK_WEBHOOK
    else
      curl -X POST -H 'Content-type: application/json' --data '{"text":"No new URS to register"}' \$LITSCAN_SLACK_WEBHOOK
    fi
    """
}

workflow search_new_ids {
    emit: new_ids
    main:
      Channel.of("Starting LitScan pipeline") | slack_message

      Channel.fromPath('workflows/litscan/queries/*.sql') \
      | get_ids \
      | check_ids \
      | sort_ids \
      | prepare_to_submit \
      | collect \
      | register_ids

      register_ids.out | register_urs

      new_ids = register_ids.out.mix(register_urs.out)
}

workflow {
  search_new_ids()
}
