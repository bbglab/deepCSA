process COMPUTE_MUTATED_EPITHELIUM {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(mutations), path(readcounts)

    output:
    tuple val(meta), path("*.exon.mutated_cells_reads.tsv")      , emit: mutated_epi_exon
    tuple val(meta), path("*.gene.mutated_cells_reads.tsv")      , emit: mutated_epi_gene
    tuple val(meta), path("*.sample.mutated_cells_reads.tsv")    , emit: mutated_epi_sample
    path  "versions.yml"                                        , topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    mutgenomes_reads_compute.py \\
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
