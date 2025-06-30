process SELECT_MUTRATES {
    tag "$meta.id"

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutation_rates)

    output:
    tuple val(meta), path("*.gene_mutrates.tsv") , emit: mutrate
    path  "versions.yml"                         , topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = task.ext.mode ?: "mutations"
    """
    omega_select_mutrate.py \\
                --mutrates ${mutation_rates} \\
                --output ${prefix}.gene_mutrates.tsv \\
                --mode ${mode};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.gene_mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
