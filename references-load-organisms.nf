nextflow.enable.dsl=2

process get_organisms {
    input:
    val(organisms)

    output:
    path("organism_textmining_mentions")

    script:
    """
    wget $organisms -O organism_textmining_mentions
    """
}

process create_csv {
    publishDir "$baseDir/workflows/references/organisms/", mode: 'copy'

    input:
    file(organism_textmining_mentions)

    output:
    path("organism_pmid")

    script:
    """
    references-create-organisms-file.py $organism_textmining_mentions organism_pmid
    """
}

process import_organisms {
    input:
    path(organism_pmid)
    path(ctl)

    """
    pgloader --on-error-stop $ctl
    """
}

workflow {
    Channel.of("https://download.jensenlab.org/organism_textmining_mentions.tsv") \
    | get_organisms \
    | create_csv | set{organism_pmid}

    load_ctl = Channel.of("$baseDir/workflows/references/organisms/load-organisms.ctl")
    import_organisms(organism_pmid, load_ctl)
}
