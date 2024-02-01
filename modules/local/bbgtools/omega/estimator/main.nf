process OMEGA_ESTIMATOR {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // TODO create a container for omega, for both the preprocessing and the estimation

    input:
    tuple val(meta), path(mutations_table), path(mutabilities_table), path(depths)
    path (annotated_panel)

    output:
    tuple val(meta), path("output_*.tsv"), emit: result
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def option = task.ext.option ?: "bayes"
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    cat > input_estimator.json << EOF
    {
    "observed_mutations_file": "${mutations_table}",
    "mutability_file":         "${mutabilities_table}",
    "depths_file"       :      "${depths}",
    "vep_annotation_file":     "${annotated_panel}",
    "grouping_folder":         "groups_dir/"
    }
    EOF
    python ../src/estimator/main.py --option ${option} --cores ${task.cpus} input_estimator.json output_${option}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def option = task.ext.option ?: "bayes"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch output_${option}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

