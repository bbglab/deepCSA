process SIGNATURES_PROBABILITIES {

    tag "${meta.id}"
    label 'process_low'

    label 'deepcsa_core'


    input:
    tuple val(meta), path (signature_probabilities)

    output:
    path ("*.decomposed_probabilities.tsv") , emit: signature_probs
    path "versions.yml"                     , topic: versions


    script:
    """
    ls ${signature_probabilities} > signature_probs_files.txt
    concat_sbs_probs.py --signature-probabilities signature_probs_files.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.decomposed_probabilities.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
