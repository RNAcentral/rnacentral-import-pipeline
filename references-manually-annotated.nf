nextflow.enable.dsl=2

process get_ids {
    publishDir "$baseDir/workflows/references/manually_annotated/", mode: 'copy'

    input:
    path(query)

    output:
    path('results')

    script:
    """
    psql -t -A -f $query "$PGDATABASE" > results
    """
}

process split_by_db {
    publishDir "$baseDir/workflows/references/manually_annotated/", mode: 'copy'

    input:
    file(results)

    output:
    path('from_*')

    script:
    """
    references-manually-annotated.py $results from_*
    """
}

workflow {
    Channel.fromPath('workflows/references/manually_annotated/query.sql') | get_ids | split_by_db
}
