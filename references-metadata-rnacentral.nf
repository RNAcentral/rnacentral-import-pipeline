nextflow.enable.dsl=2

process get_urs {
    publishDir "$projectDir/workflows/references/metadata/rnacentral", mode: 'copy'

    input:
    path(database)

    output:
    path("urs_${database.baseName}")

    script:
    """
    metadata-rnacentral.py $database urs_${database.baseName}
    """
}

process get_job {
    publishDir "$projectDir/workflows/references/metadata/rnacentral", mode: 'copy'

    input:
    path(database)

    output:
    path("job_${database.baseName}")

    script:
    """
    metadata-rnacentral.py $database job_${database.baseName}
    """
}



workflow {
    channel.fromPath('workflows/references/results/*.txt') | get_urs
    channel.fromPath('workflows/references/results/*.txt') | get_job
}
