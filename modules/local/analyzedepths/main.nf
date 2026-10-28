process ANALYZE_DEPTHS_GROUPS {
    tag "groups"
    label 'cpu_low'

    label 'deepcsa_core'

    input:
    path(features_table)
    tuple val(meta)     , path(average_depth_gene_sample)
    tuple val(meta2)    , path(average_depth_sample) 
    
    output:
    path("*.plot_depth_per_group.pdf")  , emit: plots_per_gene_per_group
    path  "versions.yml"                , topic: versions

    script:

    def separator = task.ext.separator ?: "comma"
    def custom_groups = task.ext.features_groups ? "--groups \"${task.ext.features_groups}\" " : ""
    def custom_genes = task.ext.features_genes ? "--custom-genes \"${task.ext.features_genes}\" " : ""
    def unique_identifier = task.ext.unique_identifier ? "--unique-identifier ${task.ext.unique_identifier}" : ""
    """
    depths_group_comparison.py \\
                --table-filename ${features_table} \\
                --depth-gene-sample ${average_depth_gene_sample} \\
                --depth-sample ${average_depth_sample} \\
                --separator ${separator} \\
                ${unique_identifier} \\
                ${custom_groups} \\
                ${custom_genes} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch depths_group_comparison.plot_depth_per_group.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}