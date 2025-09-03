process PLOT_OMEGA {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(mutations), path(omegas)

    output:
    tuple val(meta), path("**.pdf")  , optional: true   , emit: plots
    path "versions.yml"                                 , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    mkdir ${prefix}.plots
    plot_selection_omega.py \\
                    --sample_name ${prefix} \\
                    --mut_file ${mutations} \\
                    --omega_file ${omegas} \\
                    --outdir ${prefix}.plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def output_prefix = task.ext.output_prefix ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}${output_prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
