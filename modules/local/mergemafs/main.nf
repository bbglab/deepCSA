process MERGE_BATCH {
    tag "$meta.id"

    label 'process_high_memory'
    label 'time_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mafs)

    output:
    tuple val(meta), path("*.cohort.tsv.gz") , emit: cohort_maf
    path "versions.yml"                      , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_cohort.py --output_file ${prefix}.cohort.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch all_samples.cohort.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
