process SORT_MERGED_PANEL {

    tag "${meta.id}"

    label 'mem_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(panel)

    output:
    tuple val(meta), path("*.sorted.tsv") , emit: sorted
    path "versions.yml"                   , topic: versions

    script:
    // Sort by chromosome (field 1) and position (field 2). Assumes header in first line.
    // Using version sort for chromosome (handles chr1 chr2 chr10) after stripping 'chr' if present.
    """
    echo "[SORT_MERGED_PANEL] Sorting panel for ${meta.id}"
    head -n 1 ${panel} > sorted.tmp
    tail -n +2 ${panel} | awk 'BEGIN{OFS="\\t"} {sub(/^chr/,"",\$1); print}' | sort -k1,1V -k2,2n | awk 'BEGIN{OFS="\\t"} {print "chr"\$0}' >> sorted.tmp
    mv sorted.tmp ${panel.getBaseName()}.sorted.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${panel.getBaseName()}.sorted.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """
}
