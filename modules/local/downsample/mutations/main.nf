process DOWNSAMPLE_MUTATIONS {
    // this should be left as false since it is not deterministic
    // mutations should be downsampled differently in different runs
    cache false

    tag "$meta.id"
    label 'process_high'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutations_file)

    output:
    tuple val(meta), path("*.downsampled.mutations.tsv") , emit: downsampled_muts
    path "versions.yml"                                  , topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def downsample_prop = task.ext.downsample_prop ?: 1
    """
    downsample_script.py mutations \\
                        --file ${mutations_file} \\
                        --proportion ${downsample_prop} \\
                        --samplename ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
