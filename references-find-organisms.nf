nextflow.enable.dsl=2

process get_organisms {
    publishDir "$baseDir/workflows/references/organisms/", mode: 'copy'

    input:
    val(_flag)
    path(query)

    output:
    path("organism_pmcid")

    script:
    """
    psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > organism_pmcid
    """
}

process save_organisms {
    input:
    path(organism_pmcid)
    path(ctl)

    output:
    val('done')

    """
    pgloader --on-error-stop $ctl
    """
}

workflow find_organisms {
    take: ready
    emit: done
    main:
      query = Channel.of("$baseDir/workflows/references/organisms/get-organisms.sql")
      get_organisms(ready, query) | set{ organism_pmcid }

      save_ctl = Channel.of("$baseDir/workflows/references/organisms/save-organisms.ctl")
      save_organisms(organism_pmcid, save_ctl) | set{ done }
}

workflow {
  find_organisms()
}
