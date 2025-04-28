process DOWNSAMPLE_DEPTHS {

    tag "$meta.id"
    label 'process_high'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(depths_file)


    output:
    tuple val(meta), path("*.downsampled.tsv.gz") , emit: downsampled_depths
    path "versions.yml"                           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def downsample_prop = task.ext.downsample_prop ?: 1
    """
    downsample_script.py depths \\
                        --file ${depths_file} \\
                        --proportion ${downsample_prop}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.downsampled.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
