nextflow.enable.dsl=2

process find_organisms {
    publishDir "$baseDir/workflows/references/organisms/", mode: 'copy'

    input:
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

    """
    pgloader --on-error-stop $ctl
    """
}

workflow {
    save_ctl = Channel.of("$baseDir/workflows/references/organisms/save-organisms.ctl")
    Channel.of("$baseDir/workflows/references/organisms/find-organisms.sql") | find_organisms | set{organism_pmcid}
    save_organisms(organism_pmcid, save_ctl)
}
