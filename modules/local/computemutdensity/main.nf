process MUTATION_DENSITY {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutations), path(depth)
    tuple val(meta2), path(consensus_panel)

    output:
    tuple val(meta), path("*.mutdensities.tsv"), emit: mutdensities
    path "versions.yml", topic: versions

    script:
    def sample_name = "${meta.id}"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    compute_mutdensity.py \\
                --maf_path ${mutations} \\
                --depths_path ${depth} \\
                --annot_panel_path ${consensus_panel} \\
                --sample_name ${sample_name} \\
                --panel_version ${panel_version};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.${panel_version}.mutdensities.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
