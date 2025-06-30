process MUTATIONS_2_SIGNATURES {

    tag "${meta.id}"

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"


    input:
    tuple val(meta), path (maf), path (signature_probabilities)

    output:
    tuple val(meta), path ("*.sigs.annotated.tsv.gz")   , emit: mafs_sigs_info
    path "versions.yml"                                 , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sigs.annotated.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
