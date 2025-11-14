process COMPUTE_MATRIX {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(mut_files)

    output:
    tuple val(meta), path("*.matrix.tsv")                             , emit: matrix
    tuple val(meta), path("*.single.sigprofiler")      , optional:true, emit: single_sigprof
    tuple val(meta), path("*.per_sample")              , optional:true, emit: per_sample
    tuple val(meta), path("*.per_sample.sigprofiler")  , optional:true, emit: per_sample_sigprof
    path "versions.yml"                                               , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    mut_profile.py matrix \\
                    --sample_name ${prefix} \\
                    --mut_file ${mut_files} \\
                    --out_matrix ${prefix}.matrix.tsv \\
                    ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
