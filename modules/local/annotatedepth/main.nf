process ANNOTATE_DEPTHS {

    tag "${meta.id}"
    label 'process_low'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta) , path(depths)
    tuple val(meta2), path(panel_all)
    path (json_groups)

    output:
    // tuple val(meta), path("*.depths.annotated.tsv.gz") , emit: annotated_depths
    path("*.depths.annotated.tsv.gz")                           , emit: annotated_depths
    tuple val(meta), path("all_samples_indv.depths.tsv.gz")     , emit: all_samples_depths
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO check if the file is compressed and uncompress if needed before subsetting
    """
    cut -f 1,2,9 ${panel_all} | uniq > ${panel_all}.contexts
    merge_annotation_depths.py \\
        --annotation ${panel_all}.contexts \\
        --depths ${depths} \\
        --json_file ${json_groups}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch all_samples.cohort.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
