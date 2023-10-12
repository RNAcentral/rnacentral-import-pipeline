nextflow.enable.dsl=2

process get_expert_db_articles {
    input:
    val(_flag)
    path(query)

    output:
    path("results")

    script:
    """
    psql -t -A -f $query "$PGDB_EMBASSY_USER" > results
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
    publishDir "$baseDir/workflows/litscan/manually_annotated/", mode: 'copy'

    input:
    file(sorted_results)

    output:
    path("manually_annotated_articles")

    script:
    """
    litscan-manually-annotated.py "$PSYCOPG_CONN" $sorted_results manually_annotated_articles
    """
}

process import_manually_annotated_articles {
    input:
    path(manually_annotated_articles)
    path(ctl)

    output:
    val('done')

    """
    pgloader --on-error-stop $ctl
    """
}

workflow find_manually_annotated {
    take: ready
    emit: done
    main:
      query = Channel.fromPath('workflows/litscan/manually_annotated/query.sql')
      get_expert_db_articles(ready, query) \
      | sort_expert_db_articles \
      | find_manually_annotated_articles \
      | set{ manually_annotated_articles }

      load_ctl = Channel.of("$baseDir/workflows/litscan/manually_annotated/save-manually-annotated.ctl")
      import_manually_annotated_articles(manually_annotated_articles, load_ctl) | set{ done }
}

workflow {
  find_manually_annotated()
}
