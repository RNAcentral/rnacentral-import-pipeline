nextflow.enable.dsl=2

process create_metadata {
    input:
    path(database)

    output:
    path("metadata_${database.baseName}")

    script:
    """
    metadata.py $database metadata_${database.baseName}
    """
}

process merge_metadata {
    input:
    file(results)

    output:
    path("merged_metadata")

    script:
    """
    if [[ -f $baseDir/workflows/references/metadata/rfam/extra-ids.txt ]]; then
      cat $baseDir/workflows/references/metadata/rfam/extra-ids.txt > merged_metadata
    fi

    cat $results | sort -fb | uniq -i >> merged_metadata
    """
}

process create_xml {
    publishDir "$baseDir/workflows/references/metadata/", mode: 'copy'

    input:
    file(merged_metadata)

    output:
    path("metadata_*")

    script:
    """
    create_xml_metadata.py $merged_metadata metadata_*
    """
}

workflow {
    Channel.fromPath('workflows/references/results/*.txt') | create_metadata | merge_metadata | create_xml
}
