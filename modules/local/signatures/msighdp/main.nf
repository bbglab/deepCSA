////
// THIS SCRIPT IS CURRENTLY NOT IN USE
////


process MSIGHDP {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/msighdp:latest'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("**.pdf")     , emit: plots
    tuple val(meta), path("**.csv")     , emit: stats
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def k_guess = task.ext.k_guess ?: "12"
    """
    signatures_msighdp_run.R \\
                123 \\
                ${matrix} \\
                output.${prefix} \\
                ${k_guess} \\
                ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        mSigHdp : 2.1.2
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        mSigHdp : 2.1.2
    END_VERSIONS
    """
}
