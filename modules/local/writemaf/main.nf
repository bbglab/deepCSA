process WRITE_MAFS {

    tag "${meta.id}"
    label 'process_high_memory'

    label 'deepcsa_core'

    input:
    tuple val(meta), path(maf)
    path (json_groups)

    output:
    path("*.filtered.tsv.gz")                   , emit: mafs
    path("*.flagged-pos.bed")                   , emit: sample_flagged_beds
    path("all_samples.all-flagged-pos.bed")
    path("all_samples.cohort-wide-flagged-pos.bed")
    path "versions.yml"                         , topic: versions

    script:
    def filters = task.ext.filters ? "--filters \"${task.ext.filters}\"" : ""
    def somatic_filters = task.ext.somatic_filters ? "--somatic-filters \"${task.ext.somatic_filters}\"" : ""
    """
    write_mafs.py \\
        --maf-file ${maf} \\
        --groups-json ${json_groups} \\
        ${filters} \\
        ${somatic_filters}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch all_samples.filtered.tsv.gz
    touch all_samples.all-flagged-pos.bed
    touch all_samples.cohort-wide-flagged-pos.bed
    touch *.flagged-pos.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
