nextflow.enable.dsl=2

process get_expert_db_articles {
    input:
    path(query)

    output:
    path("results")

    script:
    """
    psql -t -A -f $query "$PGDATABASE" > results
    """
}

process sort_expert_db_articles {
    input:
    file(results)

    output:
    path("sorted_results")

    script:
    """
    cat $results | sort -fb | uniq -i >> sorted_results
    """
}

process find_manually_annotated_articles {
    publishDir "$baseDir/workflows/references/manually_annotated/", mode: 'copy'

    input:
    file(sorted_results)

    output:
    path("manually_annotated_articles")

    script:
    """
    references-manually-annotated.py "$PGDB_EMBASSY_USER" $sorted_results manually_annotated_articles
    """
}

process import_manually_annotated_articles {
    input:
    path(manually_annotated_articles)
    path(ctl)

    """
    pgloader --on-error-stop $ctl
    """
}

workflow {
    Channel.fromPath('workflows/references/manually_annotated/query.sql') \
    | get_expert_db_articles \
    | sort_expert_db_articles \
    | find_manually_annotated_articles | set{manually_annotated_articles}

    load_ctl = Channel.of("$baseDir/workflows/references/manually_annotated/save-manually-annotated.ctl")
    import_manually_annotated_articles(manually_annotated_articles, load_ctl)
}
