process PLOT_MUTDENSITY_QC {

    tag "all"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    path(all_mutdensities)
    tuple val(meta2), path(panel_file)
    path (groups_json)
    val (group_name)

    output:
    path("**.pdf")  , emit: plots
    path("**.csv")  , emit: tables

    path "versions.yml"              , topic: versions


    script:

    """
    mkdir ${group_name}.plots
    mutation_densities_qc.py \\
                    --input-file ${all_mutdensities} \\
                    --output-dir ${group_name}.plots \\
                    --panel ${panel_file} \\
                    --group-definition ${groups_json} \\
                    --group-name ${group_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
