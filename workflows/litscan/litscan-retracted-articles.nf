nextflow.enable.dsl=2

process fetch_pmcoa_metadata {
    queue 'datamover'
    memory 8.GB
    container ''

    output:
    path('out')

    """
    cp /nfs/ftp/public/databases/pmc/PMCOALiteMetadata/PMCOALiteMetadata.tgz .
    tar xvf PMCOALiteMetadata.tgz
    """
}

process check_articles {
    input:
    val(_flag)
    path(xml_dir)

    output:
    val('done')

    script:
    """
    litscan-retracted-articles.py "$PSYCOPG_CONN" $LITSCAN_SLACK_WEBHOOK ${xml_dir}
    """
}

workflow find_retracted_articles {
    take: ready
    emit: done
    main:
      fetch_pmcoa_metadata() | set { xml_dir }
      check_articles(ready, xml_dir) | set { done }
}

workflow {
    find_retracted_articles()
}
