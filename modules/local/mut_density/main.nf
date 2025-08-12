process MUTATION_DENSITY {

    tag "$meta.id"
    label 'process_high'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(somatic_mutations_file), path(depths_file), path(mutability_file)
    tuple val(meta2), path(panel_file)
    path(trinucleotide_counts_file)


    output:
    tuple val(meta), path("*.mutdensities.tsv") ,    emit: mutdensities
    path("*.logfoldchangeplot.pdf") ,                emit: mutdensities_plots
    path "versions.yml" ,                            topic: versions

    script:
    def sample_name = "${meta.id}"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    mut_density.py \\
                        --sample_name ${sample_name} \\
                        --depths_file ${depths_file} \\
                        --somatic_mutations_file ${somatic_mutations_file} \\
                        --mutability_file ${mutability_file} \\
                        --panel_file ${panel_file} \\
                        --trinucleotide_counts_file ${trinucleotide_counts_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.mutdensities.tsv
    touch ${prefix}.logfoldchangeplot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}