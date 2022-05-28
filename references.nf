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
    cat $output | sort | uniq > ${database}.txt
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
    get_unique_ids.sh ${database}.txt $database
    """
}

process submit_ids {
    input:
    tuple val(database), file("${database}_ids.txt")

    script:
    """
    upload_ids.sh ${database}_ids.txt
    """
}

workflow {
    Channel.fromPath('workflows/references/queries/*.sql') | get_ids | check_ids | sort_ids | prepare_to_submit
    // Channel.fromPath('workflows/references/queries/*.sql') | get_ids | check_ids | sort_ids | prepare_to_submit | submit_ids
}
