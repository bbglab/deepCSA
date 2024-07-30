process EXP_MUTRATE {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/axelrosendahlhuber/expected_mutrate:latest'

    input:
    tuple val(meta), path(regions)
    tuple val(meta2), path(mutations)
    tuple val(meta3), path(depths)


    output:
    tuple val(meta), path("**.png")     , emit: plots
    tuple val(meta), path("**.tsv")     , emit: stats
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def k_guess = task.ext.k_guess ?: "12"
    """
    mkdir expected_mutrate
    mutrisk_deepCSA.R \\
            $regions \\
            $mutations \\
            $depths \\
            expected_mutrate \\
            /workspace/datasets/transfer/ferriol_deepcsa/Biomart_Bladder_genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
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
    END_VERSIONS
    """
}
