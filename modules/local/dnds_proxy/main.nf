process DNDS_PROXY {
    tag "$meta.id"
    label 'process_single'

    label 'deepcsa_core'

    input:
    path(all_mutation_densities)
    tuple val(meta), path(cohort_synonymous_mutdensities)

    output:
    tuple val(meta), path("*.gene_mutdensities_n_dnds.tsv") , emit: mutdensity_with_dnds
    path  "versions.yml"                                    , topic: versions



    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def mode = task.ext.mode ?: "mutations"
    """
    mut_density_adjusted_dnds.py \\
                --mutdensities ${all_mutation_densities} \\
                --cohort-syn-mutdensities ${cohort_synonymous_mutdensities} \\
                --output ${prefix}.gene_mutdensities_n_dnds.tsv \\
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
