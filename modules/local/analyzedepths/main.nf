process ANALYZE_DEPTHS_GROUPS {

    tag "groups"
    label 'process_low'

    label 'deepcsa_core'

    input:
    path(features_table)
    // note that these input depth files are generated in another step and defined in /deepCSA/subworkflows/local/plotdepths/main.nf
    tuple val(meta) , path(average_depth_gene_sample) // needs to be added as a tupple since in PLOT_DEPTHS module (/modules/plot/depths_summary/main.nf) the output is set up as a tupple to"track" to which metadata belongs to this file
    tuple val(meta) , path(average_depth_sample) 
    


    output:
    // the main outputs will be the PDFs
    path("*.plot_depth_per_group.pdf")                           , emit: plots_per_gene_per_group

    script:

    def separator = task.ext.separator ?: "comma"
    def custom_groups = task.ext.features_groups ? "--groups \"${task.ext.features_groups}\" " : ""
    def custom_genes = task.ext.features_genes ? "--custom-genes \"${task.ext.features_genes}\" " : ""
    def unique_identifier = task.ext.unique_identifier ? "--unique-identifier ${task.ext.unique_identifier}" : ""

    // depth_group_comparison.py is in bin/ and has execution permissions add shebang 
    // ${average_depth_gene_sample} comes from subworkflows/local/plotdepths/main.nf
    """

    depth_group_comparison.py \\
                --table-filename ${features_table} \\
                --depth_gene_sample ${average_depth_gene_sample} \\
                --depth_sample ${average_depth_sample} \\
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
    touch depth_group_comparison.plot_depth_per_group.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}