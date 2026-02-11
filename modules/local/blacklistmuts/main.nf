process BLACKLIST_MUTATIONS {

    tag "$meta.id"
    label 'process_low'

    label 'deepcsa_core'

    input:
    tuple val(meta), path(mutations_file)
    path blacklist_file

    output:
    tuple val(meta), path("*.bk.filtered.mutations.tsv")  , emit: mutations
    path "versions.yml"                                , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    blacklist_muts.py \\
                    --mutations_file ${mutations_file} \\
                    --blacklist_file ${blacklist_file} \\
                    --output_file ${prefix}.bk.filtered.mutations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.bk.filtered.mutations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}