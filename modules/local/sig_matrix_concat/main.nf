process MATRIX_CONCAT {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta), path(matrix_files)

    output:
    tuple val(meta), path("*.wgs.tsv")  , emit: wgs_tsv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #for file in ${matrix_files}; do
    for file in *.WGS.sigprofiler; do
        if [ -s ${prefix}.final_matrix.wgs.tsv ]; then
            paste \$file <(cut -f 2- ${prefix}.final_matrix.wgs.tsv) > final_matrix.tsv.tmp;
        else
            mv \$file final_matrix.tsv.tmp;
        fi
        mv final_matrix.tsv.tmp ${prefix}.final_matrix.wgs.tsv;
    done;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.final_matrix.wgs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
