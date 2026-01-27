process SAMPLE_FLAGGED_BED {
    tag "$meta.id"

    label 'process_low'
    label 'time_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.flagged-pos.bed")          , emit: flagged_positions
    path "versions.yml"                                 , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def filters = task.ext.filters ?: ""
    """
    cat > mutations_subset.conf << EOF
    {
        ${filters}
    }
    EOF

    sample_flagged_positions_2bed.py \\
        --maf-file ${maf} \\
        --json-filters mutations_subset.conf \\
        --sample-name ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    """
    touch ${prefix}.flagged-pos.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
