process OMEGA_V2_RUN {
    tag "omega_v2"
    label 'cpu_medium'
    label 'process_high_memory'

    label 'deepcsa_core'

    input:
    path(mutability_tables)
    path(mutations_tables)
    path(depths_tables)
    tuple val (meta2), path(consensus_panel)
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
    mkdir -p omega_v2_workspace/deepcsa_output/depths/individual
    mkdir -p omega_v2_workspace/deepcsa_output/selection/omega/preprocessing
    mkdir -p omega_v2_workspace/deepcsa_output/regions/consensuspanels

    mv ${consensus_panel} omega_v2_workspace/deepcsa_output/regions/consensuspanels/consensus.all.tsv

    for f in ${depths_tables}; do
        mv "\$f" "omega_v2_workspace/deepcsa_output/depths/individual/."
    done

    for f in ${mutability_tables}; do
        mv "\$f" "omega_v2_workspace/deepcsa_output/selection/omega/preprocessing/."
    done

    for f in ${mutations_tables}; do
        mv "\$f" "omega_v2_workspace/deepcsa_output/selection/omega/preprocessing/."
    done

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
