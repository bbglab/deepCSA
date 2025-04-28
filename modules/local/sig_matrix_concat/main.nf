process MATRIX_CONCAT {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(matrix_files)
    path (json_samples)

    output:
    path("*_matrix*.sp.tsv")    , emit: wgs_tsv
    path("*_matrix*.hdp.tsv")   , emit: wgs_tsv_hdp
    path "versions.yml"         , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls ${matrix_files} > all_files.txt;
    concat_sigprot_matrices.py \\
                --filename_of_matrices all_files.txt \\
                --samples_json_file ${json_samples} \\
                --type_of_profile ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch samples_matrix.${prefix}.sp.tsv
    touch groups_matrix.${prefix}.sp.tsv
    touch samples_matrix.${prefix}.hdp.tsv
    touch groups_matrix.${prefix}.hdp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
