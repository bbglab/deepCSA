process COMPUTE_SATURATION_KINETICS {

    tag "$meta.id"
    label 'cpu_low'

    container 'docker.io/ferriolcalvet/saturation:v0.1.0'

    input:
    tuple val(meta) , path(mutations), path(depths), path(omega_mutability), path(relative_mutability), path(relative_mutability_index)
    tuple val(meta1), path(captured_panel_rich)
    tuple val(meta2), path(expanded_panel, stageAs: "expanded_panel.tsv")
    

    output:
    tuple val(meta), path("**.tsv"), emit: table
    tuple val(meta), path("**.pdf"), emit: plots
    tuple val(meta), path("**.png"), emit: plots_png
    path "versions.yml"            , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    discovery.py \\
                    --somatic-mutations-file ${mutations} \\
                    --vep-file ${captured_panel_rich} \\
                    --consensus-panel-file ${expanded_panel} \\
                    --depths-file ${depths} \\
                    --omega-mutability-file ${omega_mutability} \\
                    --relative-mutability-file ${relative_mutability} \\
                    --group-name ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
