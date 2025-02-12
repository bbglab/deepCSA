process DOWNSAMPLE_DEPTHS {

    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta), path(depths_file)


    output:
    tuple val(meta), path("*.downsampled.tsv.gz") , emit: downsampled_depths
    path "versions.yml"                           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def downsample_prop = task.ext.downsample_prop ?: 1
    // def configuration = task.ext.use_hotspot_bed ? "${bedfile} ${expansion} 1" : 'None 0 0'
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
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
