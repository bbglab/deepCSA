process SITE_COMPARISON {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta) , path(mutations), path(mutabilities_per_site)
    tuple val(meta2), path(annotated_panel_richer)

    output:
    tuple val(meta), path("*.comparison.tsv.gz") , emit: comparisons
    path "versions.yml"                          , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def size = task.ext.size ?: "all" // other options are 'site', 'aa_change', 'aa', '3aa', '3aa_rolling' // think if is worth having 'Naa', 'Naa_rolling'
    """
    omega_comparison_per_site.py --mutations-file ${mutations} \\
                                    --panel-file ${annotated_panel_richer} \\
                                    --mutabilities-file ${mutabilities_per_site} \\
                                    --size ${size} \\
                                    --output-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch mutabilities_per_site.${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

