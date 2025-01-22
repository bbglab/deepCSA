process SITE_COMPARISON {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/omega:20250121'

    input:
    tuple val(meta) , path(mutations), path(mutabilities_table)
    tuple val(meta2), path(annotated_panel_richer)
    // tuple val(meta2), path(annotated_panel)
    // path (genes_json)

    output:
    tuple val(meta), path("site_comparison.*.tsv.gz") , emit: comparisons
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def size = task.ext.size ?: "aa_change" // other options are 'site', 'aa_change', 'aa', '3aa', '3aa_rolling' // think if is worth having 'Naa', 'Naa_rolling'
    """

    mkdir groups;

    mv ${genes_json} groups/group_genes.json

    cat > groups/group_samples.json << EOF
    {
        "${meta.id}" : ["${meta.id}"]
    }
    EOF

    omega mutabilities --mutability-file ${mutabilities_table} \\
                        --depths-file ${depths} \\
                        --vep-annotation-file ${annotated_panel} \\
                        --grouping-folder ./groups/ \\
                        --output-fn mutabilities_per_site.${prefix}.tsv.gz \\
                        --cores 4
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

