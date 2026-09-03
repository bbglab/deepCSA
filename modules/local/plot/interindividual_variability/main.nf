process PLOT_INTERINDIVIDUAL_VARIABILITY {

    tag "samples"
    label 'cpu_low'

    label 'deepcsa_core'

    input:
    path (samples_json)
    path (all_groups_json)
    tuple val(meta), path(panel_file)
    path (mutdensities_file)
    path (adjusted_mutdensities_file)

    output:
    path("**.pdf")      , emit: plots
    path "versions.yml" , topic: versions


    script:
    def prefix = task.ext.prefix ?: "samples"
    // prefix = "${meta.id}${prefix}"
    """
    mkdir ${prefix}.variability_plots
    plot_explore_variability.py \\
                    --mutdensities ${mutdensities_file} \\
                    --adjusted-mutdensities ${adjusted_mutdensities_file} \\
                    --panel-regions ${panel_file} \\
                    --outdir ${prefix}.variability_plots \\
                    --samples-json ${samples_json} \\
                    --all-groups-json ${all_groups_json} \\


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
