process WRITE_MAFS {

    tag "${meta.id}"
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(maf)
    path (json_groups)

    output:
    path("*.filtered.tsv.gz")                   , emit: mafs
    path("*.flagged-pos.bed")                   , emit: sample_flagged_beds
    path("all_samples.cohort-flagged-pos.bed")
    path "versions.yml"                         , topic: versions

    script:
    def filters = task.ext.filters ?: ""
    def somatic_filters = task.ext.somatic_filters ?: ""
    """
    write_mafs.py \\
        --maf-file ${maf} \\
        --groups-json ${json_groups} \\
        --filters "${filters}" \\
        --somatic-filters "${somatic_filters}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch all_samples.filtered.tsv.gz
    touch all_samples.cohort-flagged-pos.bed
    touch *.flagged-pos.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
