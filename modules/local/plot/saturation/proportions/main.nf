process PLOT_SATURATION_PROPORTIONS {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta) , path(mutations)
    tuple val(meta1), path(panel_file)
    tuple val(meta2), path(captured_panel_rich)
    tuple val(meta3), path(expanded_panel)

    output:
    tuple val(meta), path("**.pdf"), optional : true    ,  emit: plots
    path "versions.yml"                                 , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    mkdir ${prefix}.plots
    plot_saturation_in_genes.py \\
                    --rich-panel ${captured_panel_rich} \\
                    --expanded-panel ${expanded_panel} \\
                    --consensus-panel ${panel_file} \\
                    --maf ${mutations} \\
                    --plots-dir ${prefix}.plots
    ##                    --genes
    ##                    --grouping-modes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
