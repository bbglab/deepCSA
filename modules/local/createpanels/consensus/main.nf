process CREATECONSENSUSPANELS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pybedtools=0.9.1--py38he0f268d_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pybedtools:0.9.1--py38he0f268d_0' :
            'biocontainers/pybedtools:0.9.1--py38he0f268d_0' }"

    input:
    tuple val(meta) , path(compact_captured_panel_annotation)
    tuple val(meta2), path(depths)
    val(consensus_min_depth)

    output:
    tuple val(meta), path("consensus*.tsv")         , emit: consensus_panel
    tuple val(meta), path("consensus*.bed")         , emit: consensus_panel_bed
    tuple val(meta), path("failing_consensus*.tsv") , emit: failing_consensus_panel
    path "versions.yml"                             , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    create_consensus_panel.py \\
                    ${compact_captured_panel_annotation} \\
                    all_samples.depths.tsv.gz \\
                    ${prefix} \\
                    $consensus_min_depth;

    bedtools merge \\
            -i <(
            tail -n +2 consensus.${prefix}.tsv | \\
            awk -F'\\t' '{print \$1, \$2-1, \$2}' OFS='\\t' | uniq
            ) > consensus.${prefix}.bed;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch consensus.${prefix}.tsv
    touch consensus.${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
