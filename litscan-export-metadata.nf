nextflow.enable.dsl=2

process create_metadata {
    input:
    val(_flag)
    path(database)

    output:
    path("metadata_${database.baseName}")

    script:
    """
    litscan-metadata.py $database metadata_${database.baseName}
    """
}

process merge_metadata {
    input:
    file(results)

    output:
    path("merged_metadata")

    script:
    """
    if [[ -f $baseDir/workflows/litscan/metadata/rfam/extra-ids.txt ]]; then
      cat $baseDir/workflows/litscan/metadata/rfam/extra-ids.txt > merged_metadata
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
    litscan-create-xml-metadata.py $merged_metadata metadata_*
    """
}

process create_release_file {
    publishDir "$params.litscan_index", mode: 'copy'

    input:
    val(_flag)

    output:
    path("release_note.txt")

    script:
    """
    litscan-create-release-note-file.sh $params.litscan_index $params.release_version
    """
}

workflow export_metadata {
    take: ready
    main:
      database = Channel.fromPath('workflows/litscan/results/*.txt')
      create_metadata(ready, database)
      | collect
      | merge_metadata
      | create_xml
      | create_release_file
}

workflow {
  export_metadata()
}
