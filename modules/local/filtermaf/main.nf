process FILTER_BATCH {
    tag "$meta.id"

    label 'process_high_memory'
    label 'time_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.cohort.filtered.tsv.gz") , emit: cohort_maf
    path "versions.yml"                               , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def repetitive_variant = task.ext.repetitive_variant ?: "${params.repetitive_variant_thres}"
    def germline_threshold = task.ext.germline_threshold ?: "${params.germline_threshold}"
    def proportion_samples_nrich = task.ext.prop_samples_nrich ?: "${params.prop_samples_nrich}"
    """
    filter_cohort.py ${maf} ${prefix} ${repetitive_variant} ${germline_threshold} ${proportion_samples_nrich}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch all_samples.cohort.filtered.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
