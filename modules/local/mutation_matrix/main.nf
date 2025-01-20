process COMPUTE_MATRIX {

    tag "$meta.id"
    label 'process_low'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta), path(mut_files)

    output:
    tuple val(meta), path("*.matrix.tsv")                             , emit: matrix
    tuple val(meta), path("*.single.sigprofiler")      , optional:true, emit: single_sigprof
    tuple val(meta), path("*.per_sample")              , optional:true, emit: per_sample
    tuple val(meta), path("*.per_sample.sigprofiler")  , optional:true, emit: per_sample_sigprof
    path "versions.yml"                                               , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
