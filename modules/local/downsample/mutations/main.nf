process DOWNSAMPLE_MUTATIONS {

    tag "$meta.id"
    label 'process_high'

    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta), path(mutations_file)

    output:
    tuple val(meta), path("*.downsampled.mutations.tsv") , emit: downsampled_muts
    path "versions.yml"                                  , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def downsample_prop = task.ext.downsample_prop ?: 0.5
    // def configuration = task.ext.use_hotspot_bed ? "${bedfile} ${expansion} 1" : 'None 0 0'
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
