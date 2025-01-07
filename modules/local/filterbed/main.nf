process FILTERBED {
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
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta)     , path(maf)
    tuple val(meta2)    , path(bedfile)

    output:
    tuple val(meta), path("*.tsv.gz")  , emit: maf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filtername = task.ext.filtername ?: "covered"
    """
    filterbed.py ${maf} ${bedfile} ${filtername};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
