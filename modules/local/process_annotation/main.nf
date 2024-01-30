process POSTPROCESS_VEP_ANNOTATION {

    tag "${vep_annotated_file}"

    label 'cpu_single'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/ferriolcalvet/bgreference"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //         'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //         'biocontainers/pandas:1.5.2' }"


    input:
    tuple val(meta), path(vep_annotated_file)

    output:
    tuple val(meta), path ("*.compact.tsv") , emit: compact_panel_annotation
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO
    // change panel postprocessing annotation into the same post processing annotation as before
    """
    zegrep -v '^##' ${vep_annotated_file} | gzip > ${vep_annotated_file}.tmp.gz
    panel_postprocessing_annotation.py \\
                    ${vep_annotated_file}.tmp.gz \\
                    ${vep_annotated_file.getBaseName()}.compact.tsv;
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${vep_annotated_file.getBaseName()}.compact.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
