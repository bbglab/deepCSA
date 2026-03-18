process MUTATIONS_2_SIGNATURES {

    tag "${meta.id}"

    label 'deepcsa_core'


    input:
    tuple val(meta), path (maf), path (signature_probabilities)

    output:
    tuple val(meta), path ("*.sigs.annotated.tsv.gz")   , emit: mafs_sigs_info
    path "versions.yml"                                 , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    signatures_mutations_n_sbs.py --mutations ${maf} \\
                            --signature-probabilities ${signature_probabilities} \\
                            --output ${prefix}.sigs.annotated.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.sigs.annotated.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
