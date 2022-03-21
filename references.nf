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
    tuple val(database), path("${database}_sorted")

    script:
    """
    cat $output | sort | uniq > ${database}_sorted
    """
}

process submit_ids {
    input:
    tuple val(database), file("${database}_sorted")

    script:
    """
    upload_ids.sh ${database}_sorted $database
    """
}

workflow {
    Channel.fromPath('workflows/references/queries/*.sql') | get_ids | check_ids | sort_ids
    // Channel.fromPath('workflows/references/queries/*.sql') | get_ids | check_ids | sort_ids | submit_ids
}
