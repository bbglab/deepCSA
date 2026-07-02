process PLOT_OMEGA_VS_DNDSCV {

    tag "${group_name}"
    label 'process_low'

    label 'deepcsa_core'

    input:
    path (all_omegas)
    path (groups_json)
    val (group_name)
    path (compiled_flagged)

    output:
    path("**.pdf")                              , optional: true    , emit: plots
    path "versions.yml" , topic: versions

    script:
    """
    mkdir ${group_name}.plots
    omega_vs_dndscv_qc.py \\
                    --input-omega-file ${all_mutdensities} \\
                    --input-dndscv-file ${} \\
                    --output-dir ${group_name}.plots \\
                    --flagged-genes-omega ${compiled_flagged}


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
