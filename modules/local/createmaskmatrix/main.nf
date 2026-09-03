process CREATE_MASK_MATRIX {
    tag "mask_matrix"

    label 'cpu_low'
    label 'mem_low'

    container "docker.io/bbglab/deepcsa-core:0.1.0"

    input:
    path(bed_files)  // List of sample-specific flagged position BED files

    output:
    path("flagged_positions.mask.tsv.gz")          , emit: mask_matrix
    path "versions.yml"                            , topic: versions

    script:
    def bed_args = bed_files.collect { "\"${it}\"" }.join(' --bed-files ')
    """
    create_mask_matrix.py \\
        --bed-files ${bed_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch flagged_positions.mask.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
