nextflow.enable.dsl=2

process get_ids {
    input:
    path(database)

    output:
    tuple val("$database.baseName"), path('results')

    script:
    """
    psql -t -A -f $database "$PGDATABASE" > results
    """
}

process check_ids {
    input:
    tuple val(database), file(results)

    output:
    tuple val(database), path('output')

    script:
    """
    check_ids.py $database $results output
    """
}

process sort_ids {
    publishDir "$baseDir/workflows/references/results/", mode: 'copy'

    input:
    tuple val(database), file(output)

    output:
    tuple val(database), path("${database}.txt")

    script:
    """
    cat $output | sort -fb | uniq -i > ${database}.txt
    """
}

process prepare_to_submit {
    publishDir "$baseDir/workflows/references/submit/", mode: 'copy'

    input:
    tuple val(database), path("${database}.txt")

    output:
    tuple val(database), path("${database}_ids.txt")

    script:
    """
    # make a copy of the old version before creating the new file
    rm -f $baseDir/workflows/references/submit/previous-release/${database}_ids.txt
    mv $baseDir/workflows/references/submit/${database}_ids.txt $baseDir/workflows/references/submit/previous-release
    get_unique_ids.sh ${database}.txt $database
    """
}

process submit_ids {
    input:
    tuple val(database), file("${database}_ids.txt")

    script:
    """
    # submit new ids only
    comm -13 $baseDir/workflows/references/submit/previous-release/${database}_ids.txt $baseDir/workflows/references/submit/${database}_ids.txt > new_${database}_ids.txt
    upload_ids.sh new_${database}_ids.txt
    """
}

workflow {
    Channel.fromPath('workflows/references/queries/*.sql') | get_ids | check_ids | sort_ids | prepare_to_submit
    // Channel.fromPath('workflows/references/queries/*.sql') | get_ids | check_ids | sort_ids | prepare_to_submit | submit_ids
}
