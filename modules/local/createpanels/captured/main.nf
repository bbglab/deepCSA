process CREATECAPTUREDPANELS {
    tag "$meta.id"
    label 'process_single'
    label 'process_medium_high_memory'

    conda "bioconda::pybedtools=0.9.1--py38he0f268d_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pybedtools:0.9.1--py38he0f268d_0' :
            'biocontainers/pybedtools:0.9.1--py38he0f268d_0' }"


    input:
    tuple val(meta), path(compact_captured_panel_annotation)

    output:
    path("*.all.tsv")                      , emit: captured_panel_all
    path("*.protein_affecting.tsv")        , emit: captured_panel_protein_affecting
    path("*.non_protein_affecting.tsv")    , emit: captured_panel_non_protein_affecting
    path("*.exons_splice_sites.tsv")       , emit: captured_panel_exons_splice_sites
    path("*.introns_intergenic.tsv")       , emit: captured_panel_introns_intergenic
    path("*.synonymous.tsv")                , emit: captured_panel_synonymous
    path("*.all.bed")                      , emit: captured_panel_all_bed
    path("*.protein_affecting.bed")        , emit: captured_panel_protein_affecting_bed
    path("*.non_protein_affecting.bed")    , emit: captured_panel_non_protein_affecting_bed
    path("*.exons_splice_sites.bed")       , emit: captured_panel_exons_splice_sites_bed
    path("*.introns_intergenic.bed")       , emit: captured_panel_introns_intergenic_bed
    path("*.synonymous.bed")                , emit: captured_panel_synonymous_bed
    path "versions.yml"                    , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_panel_versions.py \\
                    ${compact_captured_panel_annotation} \\
                    ${prefix};
    for captured_panel in \$(ls -l *.tsv | grep -v '^l' | awk '{print \$NF}'); do
        bedtools merge \\
            -i <(
                tail -n +2 \$captured_panel | \\
                awk -F'\\t' '{print \$1, \$2-1, \$2}' OFS='\\t' | uniq
            ) > \${captured_panel%.tsv}.bed;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.all.tsv
    touch ${prefix}.protein_affecting.tsv
    touch ${prefix}.non_protein_affecting.tsv
    touch ${prefix}.exons_splice_sites.tsv
    touch ${prefix}.introns_intergenic.tsv
    touch ${prefix}.all.bed
    touch ${prefix}.protein_affecting.bed
    touch ${prefix}.non_protein_affecting.bed
    touch ${prefix}.exons_splice_sites.bed
    touch ${prefix}.introns_intergenic.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
