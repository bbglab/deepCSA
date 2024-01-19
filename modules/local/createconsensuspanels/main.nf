process CREATECONSENSUSPANELS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pybedtools=0.9.1--py38he0f268d_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pybedtools:0.9.1--py38he0f268d_0' :
            'biocontainers/pybedtools:0.9.1--py38he0f268d_0' }"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //     'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(compact_captured_panel_annotation)
    tuple val(meta2), path(depths)
    val(consensus_min_depth)

    output:
    tuple val(meta), path("*.tsv"), emit: consensus_panel
    tuple val(meta), path("*.bed"), emit: consensus_panel_bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    create_consensus_panel.py \\
                    ${prefix}.compact.*.tsv \\
                    all_samples.depths.tsv.gz \\
                    $consensus_min_depth;
    for consensus_panel in \$(ls *.tsv | grep -v '^captured_panel');
    do bedtools merge \\
    -i <(tail -n +2 \$consensus_panel | \\
    awk -F'\t' '{print \$1, \$2, \$2}' OFS='\t') > \${consensus_panel%.tsv}.bed;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "TargetRegions"
    """
    touch consensus_panel.*.tsv
    touch consensus_panel.*.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
