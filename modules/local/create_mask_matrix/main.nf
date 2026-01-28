process CREATE_MASK_MATRIX {
    tag "mask_matrix"

    label 'process_low'
    label 'time_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(bed_files)  // List of sample-specific flagged position BED files

    output:
    path("flagged_positions.mask.tsv.gz")          , emit: mask_matrix
    path "versions.yml"                            , topic: versions

    script:
    """
    create_mask_matrix.py \\
        --bed-files ${bed_files.join(' --bed-files ')} \\

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
