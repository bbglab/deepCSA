process SELECT_MUTDENSITIES {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(mutation_densities)

    output:
    tuple val(meta), path("*.gene_mutdensities.tsv") , emit: mutdensity
    path  "versions.yml"                         , topic: versions



    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def mode = task.ext.mode ?: "mutations"
    """
    omega_select_mutdensity.py \\
                --mutdensities ${mutation_densities} \\
                --output ${prefix}.gene_mutdensities.tsv \\
                --mode ${mode};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.gene_mutdensities.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
