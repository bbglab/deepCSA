process FILTERBED {
    // TODO
    // reimplement it in a way that uses a BED file to know which mutations are inside
    // the regions and which ones are outside this way we avoid having to load too many
    // things in python that might slow things down
    // Look at the low mappability or low complexity filtering of the deepUMIcaller pipeline

    tag "$meta.id"

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta)     , path(maf)
    tuple val(meta2)    , path(bedfile)

    output:
    tuple val(meta), path("*.tsv.gz")  , emit: maf
    path "versions.yml"                , topic: versions

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
