process PLOTOMEGAVSOMEGAGLOBAL {

    tag "all_samples"
    label 'process_low'
    label 'deepcsa_core'

    input:
    path (all_omegas)
    path (all_omegas_globalloc)
    path (compiled_flagged)

    output:
    path("**.pdf")                              , optional: true    , emit: plots
    path("**.tsv")                              , optional: true    , emit: tables
    path "versions.yml" , topic: versions

    script:
    """
    omega_vs_omegaglobal_qc.py \\
                    --input-omega-file ${all_omegas} \\
                    --input-omegaglobal-file ${all_omegas_globalloc} \\
                    --output-dir . \\
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
