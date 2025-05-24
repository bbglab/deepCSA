process WRITE_MAFS {

    tag "${meta.id}"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(maf)
    path (json_groups)

    output:
    path("*.filtered.tsv.gz") , emit: mafs
    path "versions.yml"       , topic: versions


    script:
    // TODO reimplement with click
    """
    write_mafs.py ${maf} ${json_groups}

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
