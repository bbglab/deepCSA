process SPLIT_TSV_BY_CHROM {
    tag "$meta.id"
    label 'process_low'

    conda 'bioconda::coreutils=9.1'
    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv_chunks

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '{print > "${prefix}_"\$1".tsv"}' ${tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chr1.tsv
    touch ${prefix}_chr2.tsv
    touch ${prefix}_chr3.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1)
    END_VERSIONS
    """

}