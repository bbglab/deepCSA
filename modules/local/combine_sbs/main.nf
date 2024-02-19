process SIGNATURES_PROBABILITIES {

    tag "${meta.id}"
    label 'process_low'

    // conda "pandas:1.5.2"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //     'biocontainers/pandas:1.5.2' }"

    container 'docker.io/ferriolcalvet/bgreference'


    input:
    tuple val(meta), path (signature_probabilities)

    output:
    path ("*.decomposed_probabilities.tsv") , emit: signature_probs
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls ${signature_probabilities} > signature_probs_files.txt
    concat_sbs_probs.py --signature-probabilities signature_probs_files.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.decomposed_probabilities.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
