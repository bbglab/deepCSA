process MERGE_BATCH {
    // TODO
    // reimplement it in a way that uses a BED file to know which mutations are inside
    // the regions and which ones are outside this way we avoid having to load too many
    // things in python that might slow things down
    // Look at the low mappability or low complexity filtering of the deepUMIcaller pipeline

    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/bbglab/bgreference'

    input:
    tuple val(meta), path(mafs)

    output:
    tuple val(meta), path("*.cohort.tsv.gz") , emit: cohort_maf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_cohort.py ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch all_samples.cohort.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
