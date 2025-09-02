process PLOT_INTERINDIVIDUAL_VARIABILITY {

    tag "samples"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    path(samples_json)
    tuple val(meta), path(panel_file)
    path(mutdensities_file)

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
                    --panel-regions ${panel_file} \\
                    --outdir ${prefix}.variability_plots \\
                    --samples_json ${samples_json}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
