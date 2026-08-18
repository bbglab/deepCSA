process OMEGA_V2_RUN {
    tag "omega_v2"
    label 'cpu_medium'
    label 'process_high_memory'

    label 'deepcsa'

    input:
    path(omega_assets)
    path(context_counts)
    path(covariates)
    path(groups)


    output:
    path("data.tsv")          , emit: data
    path("omega.tsv")         , emit: omega, optional: true
    path("omega.grouped.tsv") , emit: omega_grouped
    path "versions.yml"       , topic: versions

    script:
    def args = task.ext.args ?: ""
    """
    cat > omega_v2_config.json << EOF
    {
      "path": {
        "deepcsa_output_path": "omega_v2_workspace/deepcsa_output",
        "context_counts_path": "${context_counts}",
        "covariates_path": "${covariates}",
        "samples_path": "samples.json"
      }
    }
    EOF

    run_omega_covariates.py \\
        --config omega_v2_config.json \\
        --samples-file samples.json \\
        --outfolder . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch data.tsv
    touch omega.tsv
    touch omega.grouped.tsv

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
