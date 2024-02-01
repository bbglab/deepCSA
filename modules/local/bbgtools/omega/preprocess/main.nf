process OMEGA_PREPROCESS {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // TODO create a container for omega, for both the preprocessing and the estimation

    input:
    tuple val(meta), path(mutations), path(depths), path(mutation_profile)
    path (annotated_panel)
    // TODO see if we provide directly the panel annotated as omega requires
    // see regions_sites.annotation_summary.tsv in some omega test directory

    output:
    tuple val(meta), path("mutability_per_sample_gene_context.tsv"), path("mutations_per_sample_gene_impact_context.tsv") , emit: mutabs_n_mutations_tsv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    cat > input_preprocessing.json << EOF
    {
    "vep_output_file"       : "${annotated_panel}",
    "depths_file"           : "${depths}",
    "mutations_file"        : "${mutations}",
    "mutabilities_table"    : "mutability_per_sample_gene_context.tsv",
    "table_observed_muts"   : "mutations_per_sample_gene_impact_context.tsv",
    "additional_params"     : {"mutational_profile": "${mutation_profile}"}
    }
    EOF
    python ../src/preprocessing/main.py input_preprocessing.json
    # $args -c $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 1.0
    END_VERSIONS
    """
}

