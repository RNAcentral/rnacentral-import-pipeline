nextflow.enable.dsl=2

process get_organisms {
    input:
    val(_flag)

    output:
    path("organism_textmining_mentions")

    script:
    """
    wget https://download.jensenlab.org/organism_textmining_mentions.tsv -O organism_textmining_mentions
    """
}

process create_csv {
    memory '2GB'
    publishDir "$baseDir/workflows/litscan/organisms/", mode: 'copy'

    input:
    file(organism_textmining_mentions)

    output:
    path("organism_pmid")

    script:
    """
    litscan-create-organisms-file.py $organism_textmining_mentions organism_pmid
    """
}

process import_organisms {
    input:
    path(organism_pmid)
    path(ctl)

    output:
    val('done')

    """
    pgloader --on-error-stop $ctl
    """
}

workflow load_organisms {
    take: ready
    emit: done
    main:
      get_organisms(ready) \
      | create_csv
      | set{ organism_pmid }

      load_ctl = Channel.of("$baseDir/workflows/litscan/organisms/load-organisms.ctl")
      import_organisms(organism_pmid, load_ctl) \
      | set{ done }
}

workflow {
  load_organisms()
}
