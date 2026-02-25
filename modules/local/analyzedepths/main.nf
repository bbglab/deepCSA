process ANALYZE_DEPTHS_GROUPS {

    tag "groups"
    label 'process_low'

    label 'deepcsa_core'

    input:
    path(features_table)
    //input depth file
    // note that the variable for depth is defined in /deepCSA/subworkflows/local/plotdepths/main.nf
    tuple val(meta) , path(average_depth_gene_sample) // needs to be added as a tupple since in PLOT_DEPTHS module (/modules/plot/depths_summary/main.nf) the output is set up as a tupple to"track" to which metadata belongs to this file
    


    output:
    // the main outputs will be the PDFs
    path("*.plot_depth_per_group.pdf")                           , emit: plots

    script:
    // Use meta.id to ensure each sample gets a unique folder/file name
    def output_path = task.workDir 
    def separator = task.ext.separator ? " --separator \"${task.ext.separator}\" " : ""
    def custom_groups = task.ext.features_groups ? "--groups \"${task.ext.features_groups}\" " : ""
    def custom_genes = task.ext.features_genes ? "--custom-genes \"${task.ext.features_genes}\" " : ""
    def unique_identifier = task.ext.unique_identifier ? "--unique-identifier ${task.ext.unique_identifier}" : ""

    // depth_group_comparison.py is in bin/ and has execution permissions add shebang 
    """

    depth_group_comparison.py \\
                --table-filename $features_table \\
                --depth-table $average_depth_gene_sample \\
                $separator \\
                $unique_identifier \\
                $custom_groups \\
                $custom_genes \\
                --output_prefix ${output_path}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch groups.json all_groups.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}