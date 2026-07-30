process PLOT_OMEGA_VS_GLOBAL_VS_DNDSCV {

    tag "all_samples"
    label 'process_low'
    label 'deepcsa_core'

    input:
    path (all_omegas)
    path (all_omegas_globalloc)
    path (dndscv_cv)
    path (compiled_flagged)
    path (groups_json)

    output:
    path("**.pdf")     , optional: true    , emit: plots
    path("**.tsv")     , optional: true    , emit: tables
    path "versions.yml" , topic: versions

    script:
    def dndscv_flag = params.dnds ? "--input-dndscv-file ${dndscv_cv}" : ""
    """
    omega_vs_global_vs_dndscv_qc.py \\
                    --input-omega-file ${all_omegas} \\
                    --input-omegaglobal-file ${all_omegas_globalloc} \\
                    ${dndscv_flag} \\
                    --output-dir . \\
                    --flagged-genes-omega ${compiled_flagged} \\
                    --defined-groups ${groups_json}


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
