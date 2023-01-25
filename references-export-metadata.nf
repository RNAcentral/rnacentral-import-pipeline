nextflow.enable.dsl=2

process create_metadata {
    input:
    val(_flag)
    path(database)

    output:
    path("metadata_${database.baseName}")

    script:
    """
    references-metadata.py $database metadata_${database.baseName}
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
    publishDir "$params.litscan_index", mode: 'copy'

    input:
    file(merged_metadata)

    output:
    path("metadata_*")

    script:
    """
    rm -fr "$params.litscan_index/metadata*"
    references-create-xml-metadata.py $merged_metadata metadata_*
    """
}

workflow export_metadata {
    take: ready
    main:
      database = Channel.fromPath('workflows/references/results/*.txt')
      create_metadata(ready, database)
      | collect
      | merge_metadata
      | create_xml
}

workflow {
  export_metadata()
}
