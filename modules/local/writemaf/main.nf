process WRITE_MAFS {

    tag "${meta.id}"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(maf)
    path (json_groups)

    output:
    path("*.filtered.tsv.gz") , emit: mafs
    path "versions.yml"       , topic: versions

    script:
    """
    write_mafs.py \\
        --maf-file ${maf} \\
        --groups-json ${json_groups}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch all_samples.cohort.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
