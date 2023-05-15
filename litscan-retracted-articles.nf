nextflow.enable.dsl=2

process check_articles {
    input:
    val(_flag)

    output:
    val('done')

    script:
    """
    litscan-retracted-articles.py "$PGDB_EMBASSY_USER" $LITSCAN_SLACK_WEBHOOK
    """
}

workflow find_retracted_articles {
    take: ready
    emit: done
    main:
      check_articles(ready) | set{ done }
}

workflow {
    find_retracted_articles()
}
