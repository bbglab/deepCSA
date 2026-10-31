process EFFECTIVE_DEPTH {

    tag "$meta.id"
    label 'deepcsa_core'
    label 'mem_low'

    input:
    tuple val(meta) , path(mutations_file)
    tuple val(meta2), path (all_samples_mutations)


    output:
    tuple val(meta), path("*.somatic.mutations.tsv"), emit: mutations
    path "versions.yml", topic: versions

    script:
    def prefix = "${meta.id}"
    def percentile = task.ext.percentile ?: 0.995
    """
    add_effective_depth.py \\
        --all-samples-file ${all_samples_mutations} \\
        --mutations-file ${mutations_file} \\
        --output-file ${prefix}.somatic.mutations.tsv \\
        --percentile ${percentile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.somatic.mutations.tsv"
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
