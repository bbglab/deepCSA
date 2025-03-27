process COMPUTE_MUTATED_EPITHELIUM {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutations), path(readcounts)

    output:
    tuple val(meta), path("*.exon.mutated_epithelium.tsv")      , emit: mutated_epi_exon
    tuple val(meta), path("*.gene.mutated_epithelium.tsv")      , emit: mutated_epi_gene
    tuple val(meta), path("*.sample.mutated_epithelium.tsv")    , emit: mutated_epi_sample
    path  "versions.yml"                                        , topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    compute_mutated_epithelium.py \\
                ${prefix} \\
                ${mutations} \\
                ${readcounts} ;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta.id}"
    """
    touch ${prefix}.${panel_version}.mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
