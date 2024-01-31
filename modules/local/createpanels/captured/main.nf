process CREATECAPTUREDPANELS {
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

    output:
    tuple val(meta), path("*.compact.all.tsv")                      , emit: captured_panel_all
    tuple val(meta), path("*.compact.protein_affecting.tsv")        , emit: captured_panel_protein_affecting
    tuple val(meta), path("*.compact.non_protein_affecting.tsv")    , emit: captured_panel_non_protein_affecting
    tuple val(meta), path("*.compact.exons_splice_sites.tsv")       , emit: captured_panel_exons_splice_sites
    tuple val(meta), path("*.compact.introns_intergenic.tsv")       , emit: captured_panel_introns_intergenic
    tuple val(meta), path("*.compact.all.bed")                      , emit: captured_panel_all_bed
    tuple val(meta), path("*.compact.protein_affecting.bed")        , emit: captured_panel_protein_affecting_bed
    tuple val(meta), path("*.compact.non_protein_affecting.bed")    , emit: captured_panel_non_protein_affecting_bed
    tuple val(meta), path("*.compact.exons_splice_sites.bed")       , emit: captured_panel_exons_splice_sites_bed
    tuple val(meta), path("*.compact.introns_intergenic.bed")       , emit: captured_panel_introns_intergenic_bed
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_panel_versions.py \\
                    ${compact_captured_panel_annotation} \\
                    ${prefix};
    for captured_panel in \$(ls *.tsv);
    do bedtools merge \\
    -i <(tail -n +2 \$captured_panel | \\
    awk -F'\\t' '{print \$1, \$2, \$2}' OFS='\\t') > \${captured_panel%.tsv}.bed;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "TargetRegions"
    """
    touch ${prefix}.compact.all.tsv
    touch ${prefix}.compact.protein_affecting.tsv
    touch ${prefix}.compact.non_protein_affecting.tsv
    touch ${prefix}.compact.exons_splice_sites.tsv
    touch ${prefix}.compact.introns_intergenic.tsv
    touch ${prefix}.compact.all.bed
    touch ${prefix}.compact.protein_affecting.bed
    touch ${prefix}.compact.non_protein_affecting.bed
    touch ${prefix}.compact.exons_splice_sites.bed
    touch ${prefix}.compact.introns_intergenic.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
