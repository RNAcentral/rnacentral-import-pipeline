nextflow.enable.dsl=2

process create_xml {
    publishDir "$baseDir/workflows/references/metadata/", mode: 'copy'

    input:
    path(database)

    output:
    path("metadata_${database.baseName}.xml.gz")

    script:
    """
    metadata.py $database metadata_${database.baseName}.xml.gz
    """
}

workflow {
    Channel.fromPath('workflows/references/results/*.txt') | create_xml
}
