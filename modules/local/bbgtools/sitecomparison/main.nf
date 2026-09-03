process SITE_COMPARISON {
    tag "$meta.id"

    label 'deepcsa_core'

    input:
    tuple val(meta) , path(mutations), path(mutabilities_per_site)
    tuple val(meta2), path(annotated_panel_richer)

    output:
    tuple val(meta), path("*.comparison.tsv.gz") , emit: comparisons
    path "versions.yml"                          , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def size = task.ext.size ?: "all" // other options are 'site', 'aa_change', 'aa', '3aa', '3aa_rolling' // think if is worth having 'Naa', 'Naa_rolling'
    def genes_subset = task.ext.genes_subset ?: ""
    target_genes = genes_subset != "" ? "--genes ${genes_subset}": ""
    """
    omega_comparison_per_site.py --mutations-file ${mutations} \\
                                    --panel-file ${annotated_panel_richer} \\
                                    --mutabilities-file ${mutabilities_per_site} \\
                                    --size ${size} \\
                                    --output-prefix ${prefix} \\
                                    ${target_genes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch mutabilities_per_site.${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

