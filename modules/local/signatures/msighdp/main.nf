process MSIGHDP {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/msighdp:latest'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("**.pdf")     , emit: plots
    tuple val(meta), path("**.csv")     , emit: stats
    path "versions.yml"                 , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def k_guess = task.ext.k_guess ?: "12"
    """
    msighdp_run.R \\
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
