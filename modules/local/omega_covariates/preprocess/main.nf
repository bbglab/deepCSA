process OMEGA_V2_PREPROCESS {
    tag "omega_v2"
    label 'process_single'
    label 'bbgregressions'

    input:
    path(mutability_tables)
    path(mutations_tables)
    path(depths_tables)
    tuple val (meta2), path(consensus_panel)

    output:
    path("omega_v2_workspace")    , emit: workspace
    path "versions.yml"           , topic: versions

    script:
    """
    mkdir -p omega_v2_workspace/deepcsa_output/depths/individual
    mkdir -p omega_v2_workspace/deepcsa_output/selection/omega/preprocessing
    mkdir -p omega_v2_workspace/deepcsa_output/regions/consensuspanels

    cp ${consensus_panel} omega_v2_workspace/deepcsa_output/regions/consensuspanels/consensus.all.tsv

    for f in ${depths_tables}; do
        sample=\$(basename "\$f" | cut -d'.' -f1)
        cp "\$f" "omega_v2_workspace/deepcsa_output/depths/individual/\${sample}.depths.annotated.tsv"
    done

    for f in ${mutability_tables}; do
        sample=\$(basename "\$f" | sed -E 's/^mutability_per_sample_gene_context\.//; s/\.tsv\$//')
        cp "\$f" "omega_v2_workspace/deepcsa_output/selection/omega/preprocessing/mutability_per_sample_gene_context.\${sample}.tsv"
    done

    for f in ${mutations_tables}; do
        sample=\$(basename "\$f" | sed -E 's/^mutations_per_sample_gene_impact_context\.//; s/\.tsv\$//')
        cp "\$f" "omega_v2_workspace/deepcsa_output/selection/omega/preprocessing/mutations_per_sample_gene_impact_context.\${sample}.tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p omega_v2_workspace/deepcsa_output/depths/individual
    mkdir -p omega_v2_workspace/deepcsa_output/selection/omega/preprocessing
    mkdir -p omega_v2_workspace/deepcsa_output/regions/consensuspanels

    touch omega_v2_workspace/deepcsa_output/regions/consensuspanels/consensus.all.tsv
    touch omega_v2_workspace/context_counts.csv
    touch omega_v2_workspace/covariates.tsv

    cat > omega_v2_config.json << EOF
    {
      "path": {
        "deepcsa_output_path": "omega_v2_workspace/deepcsa_output",
        "context_counts_path": "omega_v2_workspace/context_counts.csv",
        "covariates_path": "omega_v2_workspace/covariates.tsv"
      },
      "filter": {
        "sample": []
      }
    }
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
