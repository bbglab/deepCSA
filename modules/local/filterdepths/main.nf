process FILTERDEPTHS {
    tag "$meta.id"

    label 'deepcsa_core'

    input:
    tuple val(meta), path(depths)

    output:
    tuple val(meta), path("*.depths.tsv.gz")        , emit: depths
    path "versions.yml"                             , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"

    // Apply the same per-position mean depth filter used in COMPUTEDEPTHS.
    def filter_cmd = task.ext.minimum_depth
        ? "| awk 'NR == 1 {print; next} {sum = 0; for (i=3; i<=NF; i++) sum += \$i; mean = sum / (NF - 2); if (mean >= ${task.ext.minimum_depth}) print }'"
        : ""
    """
    zcat ${depths} \\
        ${filter_cmd} \\
        | gzip -c > ${prefix}.depths.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: "\$(awk --version 2>&1 | head -1)"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.depths.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: "\$(awk --version 2>&1 | head -1)"
    END_VERSIONS
    """
}
