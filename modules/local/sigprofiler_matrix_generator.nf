process SIGPROFILER_MATRIX_GENERATOR {

    tag "${meta.id}"
    label 'process_low'

    // container "docker.io/bbglab/sigprofilermatrixgenerator:latest"
    container "docker.io/ferriolcalvet/sigprofilermatrixgenerator:ucsd_dockerfile"

    input:
    tuple val(meta), path (vcf_files)

    output:
    path "output_*/"    , emit: sigprofiler_output
    path "versions.yml" , topic: versions


    script:
    def prefix = task.ext.prefix ?: ''
    prefix = "${meta.id}${prefix}"
    def ref_genome = task.ext.assembly ? "-r ${task.ext.assembly}" : ""
    """
    SigProfilerMatrixGenerator -i ${vcf_files} \\
                                -o output_${prefix} \\
                                ${ref_genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.decomposed_probabilities.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
