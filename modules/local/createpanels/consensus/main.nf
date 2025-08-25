process CREATECONSENSUSPANELS {
    tag "${meta.id}"
    label 'process_single'

    conda "python=3.10.17 bioconda::pybedtools=0.12.0 conda-forge::polars=1.30.0 conda-forge::click=8.2.1 conda-forge::gcc_linux-64=15.1.0 conda-forge::gxx_linux-64=15.1.0"
    container 'docker://bbglab/deepcsa_bed:latest'

    input:
    tuple val(meta), path(compact_captured_panel_annotation)
    tuple val(meta2), path(depths)

    output:
    tuple val(meta), path("consensus*.tsv"), emit: consensus_panel
    tuple val(meta), path("consensus*.bed"), emit: consensus_panel_bed
    tuple val(meta), path("failing_consensus*.tsv"), emit: failing_consensus_panel
    path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def args = task.ext.args ?: ""
    def genes_subset = task.ext.genes_subset ?: ""
    target_genes = genes_subset != "" ? "--genes ${genes_subset}" : ""
    """
    ## ${target_genes}
    create_consensus_panel.py \\
                    --compact_annot_panel_path ${compact_captured_panel_annotation} \\
                    --depths_path ${depths} \\
                    --version ${prefix} \\
                    ${args} \\
                    ${target_genes} \\
                    ;

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
