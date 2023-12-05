process ONCODRIVE3D {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // container = "${baseDir}/build/containers/oncodrive3d.sif"

    input:
    tuple val(meta), path(mutations), path(mutabilities)

    output:
    tuple val(meta), path("*.tsv")  , emit: tsv
    tuple val(meta), path("*.plots"), emit: plots
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "-s ${params.seed}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    oncodrive3D run -i ${mutations} -m ${mutabilities} -d ${params.datasets} -C ${prefix} -o ${prefix} $args -c $task.cpus

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


// process CLUSTERING {
//     debug true

//     //errorStrategy 'ignore'
//     container params.container
//     cpus params.cores
//     memory params.memory
//     maxForks params.max_running
//     publishDir params.outdir, mode:'copy'
//     tag "Clustering on $cohort"

//     input:
//     tuple val(cohort), path(inputs)

//     output:
//     path("${cohort}")
//                                                                  /// REMEMBER TO CHANGE OPTION WITH NEW IMG        S -> s      u -> c
//     script:
//     """
//     oncodrive3D run -i ${inputs[0]} -p ${inputs[1]} -d ${params.datasets} -C ${cohort} -o ${cohort} -s ${params.seed} -c ${params.cores}
//     """
// }
