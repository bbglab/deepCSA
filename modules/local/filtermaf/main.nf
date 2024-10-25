process FILTER_BATCH {
    tag "$meta.id"

    label 'process_high_memory'
    label 'time_low'

    conda "pandas:1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.cohort.filtered.tsv.gz") , emit: cohort_maf
    tuple val(meta), path("*.flagged-pos.bed")        , emit: flagged_muts
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def repetitive_variant = task.ext.repetitive_variant ?: "${params.repetitive_variant_thres}"
    def germline_threshold = task.ext.germline_threshold ?: "${params.germline_threshold}"
    """
    filter_cohort.py ${maf} ${prefix} ${repetitive_variant} ${germline_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch all_samples.cohort.filtered.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
