process PLOT_DEPTHS {

    tag "${meta.id}"
    label 'process_single'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(depth)
    tuple val(meta2), path(panel)

    output:
    tuple val(meta), path("*.pdf"), emit: plots
    tuple val(meta), path("*.avgdepth_per_sample.tsv"), emit: average_per_sample
    tuple val(meta), path("*depth*.tsv"), emit: depths
    path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    def plot_within_gene = task.ext.withingene ? "True" : "False"
    """
    plot_depths.py \\
                --sample_name ${prefix} \\
                --depth_file ${depth} \\
                --panel_bed6_file ${panel} \\
                --panel_name ${panel_version} \\
                --plot_within_gene ${plot_within_gene};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.${panel_version}.depths_info.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
